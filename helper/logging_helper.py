import inspect
import io
import logging
import logging.handlers
import multiprocessing
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime
from logging import LogRecord
from time import localtime, strftime
from typing import Generator, List, Optional, Type, Union

from rich.console import Console, ConsoleRenderable
from rich.live import Live
from rich.logging import RichHandler
from rich.text import Text
from rich.traceback import Traceback
from tqdm import tqdm
import functools


LOG_FILENAME = "routine_output.log"
LOG_LEVEL = logging.INFO
CONSOLE = Console()

GRAY = "\x1b[38;21m"
WHITE = "\x1b[38;5;15m"
YELLOW = "\x1b[38;5;226m"
RED = "\x1b[38;5;196m"
BOLD_RED = "\x1b[31;1m"
RESET = "\x1b[0m"


def is_in_ipython():
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            return True
    except ImportError:
        pass
    return False


def capture_rich_renderable_as_string(renderable, width: int = 200) -> str:
    string_io = io.StringIO()
    capture_console = Console(file=string_io, record=True, width=width)
    capture_console.print(renderable)
    return string_io.getvalue()


def worker_init(log_queue: multiprocessing.Queue, level: int = logging.INFO):
    root = logging.getLogger()
    if root.hasHandlers():
        root.handlers.clear()
    root.setLevel(level)
    handler = logging.handlers.QueueHandler(log_queue)
    root.addHandler(handler)


def grouped_logs(arg=None):
    def _apply_logging(func, extractor):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                if extractor:
                    target = extractor(*args, **kwargs)
                    if isinstance(target, logging.Logger):
                        logger = target
                        name = target.name
                    else:
                        name = str(target)
                        logger = logging.getLogger(name)
                else:
                    name = func.__name__
                    logger = logging.getLogger(name)
            except Exception:
                name = func.__name__
                logger = logging.getLogger(name)

            with LogContext(logger).grouped_logs(name):
                return func(*args, **kwargs)

        return wrapper

    if callable(arg) and getattr(arg, "__name__", "") != "<lambda>":
        return _apply_logging(arg, None)

    extractor = arg

    def decorator(func):
        return _apply_logging(func, extractor)

    return decorator


class BufferedWorkerHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)


class RealtimeInjector(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        if not getattr(record, 'summary', False):
            record.realtime = True
        return True


class ConsoleDisplayFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return not getattr(record, 'summary', False)


class FileDisplayFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return not getattr(record, 'realtime', False)


class NoFileOnlyFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        return not getattr(record, 'file_only', False)


class CustomRichHandler(RichHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def render_message(self, record: LogRecord, message: str) -> ConsoleRenderable:
        if getattr(record, 'summary', False):
            return Text.from_markup(message) if self.markup else Text(message)
        message_renderable = super().render_message(record, message)
        return Text.assemble(
            Text(f"{record.name} ", style="bold cyan"),
            message_renderable,
        )


class MarkupStrippingRichHandler(CustomRichHandler):
    def render_message(self, record: LogRecord, message: str) -> ConsoleRenderable:
        """Strips markup from the message before passing it to the parent renderer."""
        plain_message = Text.from_markup(message).plain
        return super().render_message(record, plain_message)

    def emit(self, record: LogRecord):
        if getattr(record, 'summary', False):
            try:
                msg = self.format(record)
                plain_msg = Text.from_markup(msg).plain
                self.console.file.write(plain_msg + "\n")
                self.console.file.flush()
            except Exception:
                self.handleError(record)
            return
        super().emit(record)

    def render(self, *, record: LogRecord, traceback: Optional[Traceback], message_renderable: ConsoleRenderable) -> ConsoleRenderable:
        if getattr(record, 'summary', False):
            time_format = self._log_render.time_format
            time_str = datetime.fromtimestamp(record.created).strftime(time_format)

            level_text = self.get_level_text(record)
            level_str = level_text.plain.ljust(8)
            plain_str = str(message_renderable)
            lines = plain_str.split('\n')
            indented_lines = ['        ' + line for line in lines]
            indented_str = '\n'.join(indented_lines)
            indented_message = Text(indented_str)
            return Text.assemble(Text(f"{time_str} {level_str}\n"), indented_message)
        else:
            return super().render(record=record, traceback=traceback, message_renderable=message_renderable)


class _DuplicateFilter:
    def __init__(self) -> None:
        self.msgs = set()

    def filter(self, record: logging.LogRecord) -> bool:
        if record.msg in self.msgs:
            return False
        self.msgs.add(record.msg)
        return True


def setup_logging(
    output_file: Union[str, None] = None,
    logger: logging.Logger = logging.getLogger(""),
    level: Union[int, None] = None,
    console_markup: bool = False,
    queue: Union[multiprocessing.Queue, None] = None,
) -> logging.Logger:

    if queue is not None:
        if logger.hasHandlers():
            logger.handlers.clear()
        logger.setLevel(level or LOG_LEVEL)
        handler = logging.handlers.QueueHandler(queue)
        logger.addHandler(handler)
        logger.propagate = False
        return logger

    root = logging.getLogger()
    if any(isinstance(h, logging.handlers.QueueHandler) for h in root.handlers):
        logger.setLevel(level or LOG_LEVEL)
        return logger

    if output_file is None:
        output_file = LOG_FILENAME
    if level is None:
        level = LOG_LEVEL

    if logger.hasHandlers():
        logger.handlers.clear()

    logger.setLevel(level)

    console_handler = CustomRichHandler(
        console=CONSOLE,
        rich_tracebacks=True,
        show_time=True,
        show_level=True,
        show_path=True,
        log_time_format="[%Y-%m-%d %H:%M:%S]",
        markup=console_markup,
    )
    console_handler.addFilter(NoFileOnlyFilter())
    console_handler.addFilter(ConsoleDisplayFilter())
    logger.addHandler(console_handler)

    log_file = open(output_file, "a")
    file_console = Console(file=log_file, record=True, width=172)
    file_handler = MarkupStrippingRichHandler(
        console=file_console,
        show_time=True,
        show_level=True,
        show_path=True,
        log_time_format="[%Y-%m-%d %H:%M:%S.%f]",
        rich_tracebacks=False,
        markup=False,
    )
    file_handler.addFilter(FileDisplayFilter())
    logger.addHandler(file_handler)

    # Install the duplicate filter permanently if not already present.
    if not any(isinstance(f, _DuplicateFilter) for f in logger.filters):
        logger.addFilter(_DuplicateFilter())

    return logger


class _UnclosableStream:
    def __init__(self, stream):
        self._stream = stream

    def write(self, *args, **kwargs):
        return self._stream.write(*args, **kwargs)

    def flush(self, *args, **kwargs):
        return self._stream.flush(*args, **kwargs)

    def close(self):
        pass

    def __getattr__(self, name):
        return getattr(self._stream, name)


class TqdmRichLiveIO(io.StringIO):
    def __init__(self, live_instance: Live, handler: RichHandler, logger_name: str):
        super().__init__()
        self.live = live_instance
        self.handler = handler
        self.logger_name = logger_name
        self.last_text = ""
        self._find_caller()

    def _find_caller(self):
        self.caller_filename = ""
        self.caller_lineno = 0
        for frame_info in inspect.stack():
            if 'tqdm' not in frame_info.filename and __file__ not in frame_info.filename:
                self.caller_filename = frame_info.filename
                self.caller_lineno = frame_info.lineno
                break

    def write(self, s: str):
        text = s.strip()
        if text and text != self.last_text:
            self.last_text = text
            record = logging.LogRecord(
                name=self.logger_name,
                level=logging.INFO,
                pathname=self.caller_filename,
                lineno=self.caller_lineno,
                msg=text,
                args=(),
                exc_info=None,
            )
            message_renderable = self.handler.render_message(record, text)
            full_renderable = self.handler.render(
                record=record,
                traceback=None,
                message_renderable=message_renderable
            )
            self.live.update(full_renderable)

    def flush(self):
        pass


class LogContext:

    def __init__(self, logger: logging.Logger) -> None:
        self.logger = logger

    @contextmanager
    def grouped_logs(self, worker_name: str) -> Generator[None, None, None]:
        class SimpleBuffer(logging.Handler):
            def __init__(self):
                super().__init__()
                self.buffer = []

            def emit(self, record):
                msg_text = self.format(record)

                filename = os.path.basename(record.pathname)
                path_info = f"{filename}:{record.lineno}"

                formatted_lines, lines, target_width, n_indent = [], msg_text.split("\n"), 172, 4

                for i, line in enumerate(lines):
                    indented_line = " " * n_indent + line
                    if i == 0:
                        current_len = len(indented_line)
                        needed_padding = target_width - current_len - len(path_info)
                        if needed_padding < 2:
                            needed_padding = 2
                        full_line = f"{indented_line}{' ' * needed_padding}{path_info}"
                    else:
                        full_line = indented_line
                    formatted_lines.append(full_line)

                self.buffer.append("\n".join(formatted_lines))

        buffer = SimpleBuffer()
        buffer.setFormatter(logging.Formatter(
            "[%(asctime)s.%(msecs)03d] %(levelname)-8s %(name)s: %(message)s",
            datefmt="%H:%M:%S"
        ))
        injector = RealtimeInjector()
        self.logger.addHandler(buffer)
        self.logger.addFilter(injector)
        try:
            yield
        finally:
            self.logger.removeHandler(buffer)
            self.logger.removeFilter(injector)
            if buffer.buffer:
                block = "\n".join(buffer.buffer)
                header = f"[bold cyan]{'=' * 30} START WORKER: {worker_name} {'=' * 30}[/]"
                footer = f"[bold cyan]{'=' * 31} END WORKER: {worker_name} {'=' * 31}[/]"
                self.logger.info(f"{header}\n{block}\n{footer}", extra={'summary': True})

    @contextmanager
    def parallel_session(self) -> Generator[dict, None, None]:
        manager = multiprocessing.Manager()
        log_queue = manager.Queue()

        listener = logging.handlers.QueueListener(
            log_queue,
            *self.logger.handlers,
            respect_handler_level=True
        )
        listener.start()
        pool_config = {
            "initializer": worker_init,
            "initargs": (log_queue, self.logger.level)
        }

        try:
            yield pool_config
        finally:
            listener.stop()
            manager.shutdown()

    @contextmanager
    def redirect_tqdm(self) -> Generator[None, None, None]:
        handler = None
        for h in self.logger.handlers:
            if isinstance(h, RichHandler) and hasattr(h, 'console'):
                handler = h
                break

        if not handler:
            raise ValueError("A RichHandler with a console attribute must be attached to the logger.")

        console = handler.console
        original_init = tqdm.__init__

        with Live(console=console, transient=True, refresh_per_second=20) as live:
            tqdm_io = TqdmRichLiveIO(live, handler, self.logger.name)

            def patched_init(self, *args, **kwargs):
                if 'file' not in kwargs:
                    kwargs['file'] = _UnclosableStream(tqdm_io)
                if 'disable' not in kwargs:
                    kwargs['disable'] = False

                original_init(self, *args, **kwargs)

            tqdm.__init__ = patched_init
            try:
                yield
            finally:
                tqdm.__init__ = original_init

                if tqdm_io.last_text:
                    self.logger.info(tqdm_io.last_text, stacklevel=3)

    @contextmanager
    def suppress_console_logging(self) -> Generator[None, None, None]:
        original_handlers = list(self.logger.handlers)  # Make a copy
        console_handlers_to_remove: List[logging.Handler] = []

        for handler in original_handlers:
            if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                console_handlers_to_remove.append(handler)
                self.logger.removeHandler(handler)

        try:
            yield
        finally:
            for handler in console_handlers_to_remove:
                if handler not in self.logger.handlers:  # Avoid adding duplicates if somehow re-added
                    self.logger.addHandler(handler)

            current_handlers = list(self.logger.handlers)
            for handler in original_handlers:
                if handler not in current_handlers and handler not in console_handlers_to_remove:
                    self.logger.addHandler(handler)

    @contextmanager
    def duplicate_filter(self) -> Generator[None, None, None]:
        if any(isinstance(f, _DuplicateFilter) for f in self.logger.filters):
            yield
        else:
            dup_filter = _DuplicateFilter()
            self.logger.addFilter(dup_filter)
            try:
                yield
            finally:
                self.logger.removeFilter(dup_filter)

    @contextmanager
    def logging_raised_Error(self) -> Generator[None, None, None]:
        try:
            yield
        except Exception as e:
            self.logger.error(e)
            raise

    @contextmanager
    def set_logging_level(self, level: int) -> Generator[None, None, None]:
        _old_level = self.logger.level
        self.logger.setLevel(level)
        try:
            yield
        finally:
            self.logger.setLevel(_old_level)

    @contextmanager
    def suppress_logging(self) -> Generator[None, None, None]:
        original_level = self.logger.level
        self.logger.setLevel(logging.CRITICAL + 1)
        try:
            yield
        finally:
            self.logger.setLevel(original_level)

    @contextmanager
    def log_and_suppress(self, *exceptions: Type[Exception], msg: str = "An exception was suppressed"):
        try:
            yield
        except exceptions or (Exception,) as e:
            self.logger.error(f"{msg}: {type(e).__name__} - {e}", exc_info=True)

    @contextmanager
    def suppress_terminal_print(self) -> Generator[None, None, None]:
        if is_in_ipython():
            from IPython.display import display
            from ipywidgets import Output

            out = Output()
            with out:
                yield
        else:
            with open(os.devnull, 'w') as f, redirect_stdout(f), redirect_stderr(f):
                yield
