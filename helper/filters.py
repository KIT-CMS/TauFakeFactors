from typing import Any


def emb_tau_gen_match(rdf: Any, channel: str) -> Any:
    """
    Function to make sure the tau pair final state objects originated for genuine taus for embedded events.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        rdf = rdf.Filter(
            "(gen_match_1 == 3) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(gen_match_1 == 4) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(gen_match_1 == 5) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    else:
        raise ValueError(
            f"Eventfilter: emb tau gen match: Such a channel is not defined: {channel}"
        )
    return rdf


def emb_boostedtau_gen_match(rdf: Any, channel: str) -> Any:
    """
    Function to make sure the boostedtau pair final state objects originated for genuine taus for embedded events.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 3) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 4) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    else:
        raise ValueError(
            f"Eventfilter: emb boostedtau gen match: Such a channel is not defined: {channel}"
        )
    return rdf


def tau_origin_split(rdf: Any, channel: str, tau_gen_mode: str) -> Any:
    """
    Function to apply a cut based on the origin of the tau pair. This is needed to prevent an overestimation because
    genuine tau pairs are estimated with embedded events and events with jets faking hadronic taus are the targed the
    fake factor measurement. There are diffrent possibilities where a lepton (electron, muon, hadronic tau) originate from.
    Electrons and muons are either genuine or decay from a tau. Hadtronic taus are either decays from a genuine tau
    or jets/leptons misidentified as a hadronic tau.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"
        tau_gen_mode: There are 4 options:
                      "T" both taus in the tau pair originated from genuine taus,
                      "J" a hadronic tau is faked by a jet,
                      "L" everything else where mainly the muons or electrons are genuine,
                      "all" when no split should be applied

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 3) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(gen_match_1==3 && gen_match_2==5)) && (!(gen_match_2 == 6))",
                "tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: tau gen match: et: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "mt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 4) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(gen_match_1==4 && gen_match_2==5)) && (!(gen_match_2 == 6))",
                "tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: tau gen match: mt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "tt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 5) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_1 == 6) || (gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!((gen_match_1 == 5) && (gen_match_2 == 5)) && !((gen_match_1 == 6) || (gen_match_2 == 6)))",
                "tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: tau gen match: tt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    else:
        raise ValueError(
            f"Eventfilter: tau gen match: Such a channel is not defined: {channel}"
        )

    return rdf


def boostedtau_origin_split(rdf: Any, channel: str, tau_gen_mode: str) -> Any:
    """
    Function to apply a cut based on the origin of the boosted tau pair. This is needed to prevent an overestimation because
    genuine tau pairs are estimated with embedded events and events with jets faking hadronic taus are the targed the
    fake factor measurement. There are diffrent possibilities where a lepton (electron, muon, hadronic tau) originate from.
    Electrons and muons are either genuine or decay from a tau. Hadtronic taus are either decays from a genuine tau
    or jets/leptons misidentified as a hadronic tau.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"
        tau_gen_mode: There are 4 options:
                      "T" both taus in the boosted tau pair originated from genuine taus,
                      "J" a hadronic tau is faked by a jet,
                      "L" everything else where mainly the muons or electrons are genuine,
                      "all" when no split should be applied

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 3) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(boosted_gen_match_1==3 && boosted_gen_match_2==5)) && (!(boosted_gen_match_2 == 6))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: et: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "mt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 4) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(boosted_gen_match_1==4 && boosted_gen_match_2==5)) && (!(boosted_gen_match_2 == 6))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: mt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "tt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 6) || (boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!((boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)) && !((boosted_gen_match_1 == 6) || (boosted_gen_match_2 == 6)))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: tt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    else:
        raise ValueError(
            f"Eventfilter: tau gen match: Such a channel is not defined: {channel}"
        )

    return rdf
