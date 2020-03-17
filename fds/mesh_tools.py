"""!
BlenderFDS, mesh alignment.
"""

# Mesh alignment:
#
# Before:
#   |   |rcs|   |   | ri  ref mesh
#   ·---·---·---·---·
#  rx0   <-rl->    rx1
# mx0      <-ml->      mx1
#  ·------·------·------·
#  |      | mcs  |      | mi  other mesh
# ---> axis

# Either protect rl or rcs.
# Instead ml and mcs are changed for alignment.

# After:
#  |   |   |   |   |
#  ·---·---·---·---·-------·
#  |       |       |       |

# Not allowed:
#  |    |     |    |
#  ·----·--·--·----·
#  |       |       |

#  |   |   |   |
#  ·---·---·---·---·
#  |       |       |


def _factor(n):
    """!
    Generator for prime factors of n.
    Many thanks Dhananjay Nene (http://dhananjaynene.com/)
    for publishing this code
    """
    yield 1
    i = 2
    limit = n ** 0.5
    while i <= limit:
        if n % i == 0:
            yield i
            n = n / i
            limit = n ** 0.5
        else:
            i += 1
    if n > 1:
        yield int(n)


def _n_for_poisson(n):
    """!Get a good number for poisson solver at least bigger than n."""
    good = set((1, 2, 3, 5))
    while True:
        if [i for i in _factor(n) if i not in good]:
            n += 1
        else:
            break
    return n


def _align_along_axis(ri, rx0, rx1, mi, mx0, mx1, poisson=False, protect_rl=False):
    """!Align coarse MESH to fixed ref MESH along an axis."""
    # Init
    assert rx0 < rx1 and mx0 < mx1  # coordinate order
    rl, ml = rx1 - rx0, mx1 - mx0  # lengths
    rcs, mcs = rl / ri, ml / mi  # cell sizes
    # Coarsening ratio
    assert mcs / rcs > 0.501  # same or coarser allowed, protect from float err
    n = round(mcs / rcs)
    # Set ref cell count multiple of n
    # to allow full cover of coarse cells by ref cells
    ri = round(ri / n) * n
    if poisson:
        ri = _n_for_poisson(ri)
    if protect_rl:  # protect ref length
        rcs = rl / ri  # reduce ref cell size
    else:
        rl = rcs * ri  # extend ref length due to updated ri
        rx1 = rx0 + rl
    # Calc new coarse cell size from ref cell size
    mcs = rcs * n
    # Calc new coarse cell count,
    # trying to keep ml as close as possible to the original
    mi = round(ml / mcs)
    if poisson:
        mi = _n_for_poisson(mi)
    # Align coarse mesh positions to the ref mesh
    mx0 = rx0 + round((mx0 - rx0) / mcs) * mcs
    ml = mcs * mi  # extend other mesh due to updated mi
    mx1 = mx0 + ml
    print("n:", n, "rcs:", rcs, "mcs:", mcs)  # FIXME
    return ri, rx0, rx1, mi, mx0, mx1


def _align_along_x(rijk, rxbs, mijk, mxbs, poisson=False, protect_rl=False):
    """!Align coarse MESH to fixed ref MESH along axis x."""
    rijk[0], rxbs[0], rxbs[1], mijk[0], mxbs[0], mxbs[1] = _align_along_axis(
        ri=rijk[0],
        rx0=rxbs[0],
        rx1=rxbs[1],
        mi=mijk[0],
        mx0=mxbs[0],
        mx1=mxbs[1],
        poisson=False,  # not needed along x
        protect_rl=protect_rl,
    )


def _align_along_y(rijk, rxbs, mijk, mxbs, poisson=False, protect_rl=False):
    """!Align coarse MESH to fixed ref MESH along axis y."""
    rijk[1], rxbs[2], rxbs[3], mijk[1], mxbs[2], mxbs[3] = _align_along_axis(
        ri=rijk[1],
        rx0=rxbs[2],
        rx1=rxbs[3],
        mi=mijk[1],
        mx0=mxbs[2],
        mx1=mxbs[3],
        poisson=poisson,  # needed along y
        protect_rl=protect_rl,
    )


def _align_along_z(rijk, rxbs, mijk, mxbs, poisson=False, protect_rl=False):
    """!Align coarse MESH to fixed ref MESH along axis z."""
    rijk[2], rxbs[4], rxbs[5], mijk[2], mxbs[4], mxbs[5] = _align_along_axis(
        ri=rijk[2],
        rx0=rxbs[4],
        rx1=rxbs[5],
        mi=mijk[2],
        mx0=mxbs[4],
        mx1=mxbs[5],
        poisson=poisson,  # needed along z
        protect_rl=protect_rl,
    )


# rx0 rx1
#  .---.
#       delta .----.
#            mx0  mx1


def _is_far(rxbs, mxbs, deltas):
    return (
        rxbs[0] - mxbs[1] > deltas[0]
        or mxbs[0] - rxbs[1] > deltas[0]  # x
        or rxbs[2] - mxbs[3] > deltas[1]
        or mxbs[2] - rxbs[3] > deltas[1]  # y
        or rxbs[4] - mxbs[5] > deltas[2]
        or mxbs[4] - rxbs[5] > deltas[2]  # z
    )


def align_meshes(rijk, rxbs, mijk, mxbs, poisson=False, protect_rl=False):
    """!
    Function to align meshes.
    @param rijk:  ijk of the ref mesh.
    @param rxbs: xbs of the ref mesh.
    @param mijk: ijk of the other mesh.
    @param mxbs: xbs of the other mesh.
    @param poisson: True for respecting the Poisson constraint.
    @param protect_rl: True to protect ref length.
    @return return new rijk, rxbs, mijk and mxbs.
    """
    # Init
    deltas = (  # rcs
        abs(rxbs[0] - rxbs[1]) / rijk[0],
        abs(rxbs[2] - rxbs[3]) / rijk[1],
        abs(rxbs[4] - rxbs[5]) / rijk[2],
    )
    # Are meshes far apart?
    if _is_far(rxbs=rxbs, mxbs=mxbs, deltas=deltas):
        print("Meshes are far apart, no modification")  # FIXME
        return rijk, rxbs, mijk, mxbs
    # If mesh sides are close, then snap them
    # otherwise align their meshes
    if abs(rxbs[0] - mxbs[1]) <= deltas[0]:  # -x close?
        print("Snapping -x")
        mxbs[1] = rxbs[0]
    elif abs(mxbs[0] - rxbs[1]) <= deltas[0]:  # +x close?
        print("Snapping +x")
        mxbs[0] = rxbs[1]
    else:
        print("Aligning x")
        _align_along_x(rijk, rxbs, mijk, mxbs, poisson, protect_rl)
    if abs(rxbs[2] - mxbs[3]) <= deltas[1]:  # -y close?
        print("Snapping -y")
        mxbs[3] = rxbs[2]
    elif abs(mxbs[2] - rxbs[3]) <= deltas[1]:  # +y close?
        print("Snapping +y")
        mxbs[2] = rxbs[3]
    else:
        print("Aligning y")
        _align_along_y(rijk, rxbs, mijk, mxbs, poisson, protect_rl)
    if abs(rxbs[4] - mxbs[5]) <= deltas[2]:  # -z close?
        print("Snapping -z")
        mxbs[5] = rxbs[4]
    elif abs(mxbs[4] - rxbs[5]) <= deltas[2]:  # +z close?
        print("Snapping +z")
        mxbs[4] = rxbs[5]
    else:
        print("Aligning z")
        _align_along_z(rijk, rxbs, mijk, mxbs, poisson, protect_rl)
    print("ref:", rijk, rxbs)  # FIXME
    print("oth:", mijk, mxbs)
    return rijk, rxbs, mijk, mxbs


def calc_poisson_ijk(ijk):
    """!
    Get an IJK respecting the Poisson constraint, close to the current one.
    @param ijk: ijk of the mesh.
    @return return new ijk values.
    """
    return ijk[0], _n_for_poisson(ijk[1]), _n_for_poisson(ijk[2])


def calc_cell_sizes(ijk, xbs):
    """!
    Calc MESH cell sizes.
    @param ijk: ijk of the mesh.
    @param xbs: xbs of the mesh.
    @return return the MESH cell sizes.
    """
    return (
        (xbs[1] - xbs[0]) / ijk[0],
        (xbs[3] - xbs[2]) / ijk[1],
        (xbs[5] - xbs[4]) / ijk[2],
    )


def calc_ijk(xbs, desired_cs, poisson):
    """!
    Calc MESH IJK from cell sizes.
    @param xbs: xbs of the mesh.
    @param desired_cs: desired cell sizes.
    @param poisson: True for respecting the Poisson constraint.
    @return return the MESH IJK from cell sizes.
    """
    print(xbs, desired_cs, poisson)  # FIXME
    ijk = (
        round((xbs[1] - xbs[0]) / desired_cs[0]) or 1,
        round((xbs[3] - xbs[2]) / desired_cs[1]) or 1,
        round((xbs[5] - xbs[4]) / desired_cs[2]) or 1,
    )
    if poisson:
        return calc_poisson_ijk(ijk)
    else:
        return ijk


def calc_cell_infos(ijk, xbs):
    """!
    Calc many cell infos: cell ijk goodness, sizes, count and aspect ratio.
    @param ijk: ijk of the mesh.
    @param xbs: xbs of the mesh.
    @return return if cell infos.
    """
    cs = calc_cell_sizes(ijk, xbs)
    has_good_ijk = tuple(ijk) == calc_poisson_ijk(ijk)
    cell_count = ijk[0] * ijk[1] * ijk[2]
    cell_sizes_sorted = sorted(cs)
    try:
        cell_aspect_ratio = max(
            cell_sizes_sorted[2] / cell_sizes_sorted[0],
            cell_sizes_sorted[2] / cell_sizes_sorted[1],
            cell_sizes_sorted[1] / cell_sizes_sorted[0],
        )
    except ZeroDivisionError:
        cell_aspect_ratio = 999.0
    # Return
    return has_good_ijk, cs, cell_count, cell_aspect_ratio


def test():
    print("Test")
    align_meshes(
        rijk=[15, 37, 51],
        rxbs=[0, 5, 0, 5, 0, 5],
        mijk=[9, 38, 20],
        mxbs=[0, 5, -5, 0, -5, 0],
        poisson=True,
        protect_rl=True,
    )


if __name__ == "__main__":
    test()
