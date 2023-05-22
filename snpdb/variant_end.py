
# Modified from https://github.com/jamescasbon/PyVCF/blob/master/vcf/model.py#L133

def get_start_end(position, ref, alt):

    if alt is None:
        start, end = _compute_coordinates_for_none_alt(position, ref)
    else:
        ref_length = len(ref)
        alt_length = len(alt)
        if ref_length == alt_length == 1:
            start, end = _compute_coordinates_for_snp(position, ref)
        elif "<" in alt:
            start, end = _compute_coordinates_for_sv(position, ref)
        else:
            start, end = _compute_coordinates_for_indel(position, ref)
    return start, end


def _compute_coordinates_for_none_alt(position, ref):
    start = position - 1
    end = start + len(ref)
    return (start, end)


def _compute_coordinates_for_snp(position, ref):
    if len(ref) > 1:
        start = position
        end = start + (len(ref) - 1)
    else:
        start = position - 1
        end = position
    return (start, end)


def _compute_coordinates_for_indel(position, ref):
    if len(ref) > 1:
        start = position
        end = start + (len(ref) - 1)
    else:
        start = end = position
    return (start, end)


def _compute_coordinates_for_sv(position, ref):
    start = position - 1
    end = start + len(ref)
    return (start, end)
