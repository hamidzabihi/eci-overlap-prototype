import numpy as np

def read_complex_section(path, secnum):
    """Read section '! N <title>' of SHARC QM.out as a complex numpy array."""
    lines = open(path).read().splitlines()
    start = None
    for i, L in enumerate(lines):
        if L.startswith(f'! {secnum} '):
            start = i
            break
    if start is None:
        raise ValueError(f"Section ! {secnum} not found")

    # First non-empty line after header is dimensions "r c" or "n"
    dim_line = lines[start + 1].strip()
    dims = [int(x) for x in dim_line.split()]
    if len(dims) == 1:
        shape = (dims[0],)
        n_vals = dims[0]
    else:
        shape = tuple(dims)
        n_vals = int(np.prod(dims))

    # Next n_vals * 2 floats are real, imag alternating
    vals = []
    i = start + 2
    while len(vals) < n_vals * 2 and i < len(lines):
        L = lines[i]
        if L.startswith('!'):
            break
        try:
            vals.extend(float(x) for x in L.split())
        except ValueError:
            pass
        i += 1

    arr = np.array(vals[:n_vals*2]).reshape(n_vals, 2)
    return (arr[:, 0] + 1j * arr[:, 1]).reshape(shape)


if __name__ == '__main__':
    p = '/home/hamid/sharc4/tests/interface/inputs/qmout/qmout_overlap.out'
    S = read_complex_section(p, 6)
    ph = read_complex_section(p, 7)
    print("S shape:", S.shape, "dtype:", S.dtype)
    print(S)
    print()
    print("phases shape:", ph.shape, "dtype:", ph.dtype)
    print(ph)
