function set_size(width; fraction=1.0, subplots=(1, 1))
    """
    Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predefined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: Tuple, optional
            The number of rows and columns of subplots.

    Returns
    -------
    fig_dim: Tuple
            Dimensions of figure in inches
    """
    if width == "thesis"
        width_pt = 426.79135
    elseif width == "beamer"
        width_pt = 307.28987
    else
        width_pt = width
    end

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5.0^0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[1] / subplots[2])

    return (fig_width_in, fig_height_in)
end