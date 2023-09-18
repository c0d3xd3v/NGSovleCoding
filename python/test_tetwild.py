import wildmeshing as wm

path = "/home/kai/Development/github/NGSovleCoding/data/tuningfork.stl"
tetra = wm.Tetrahedralizer(
    stop_quality=1000, 
    skip_simplify=True, 
    coarsen=False,
    edge_length_r=0.05)

tetra.load_mesh(path)
tetra.tetrahedralize()

tetra.save("/home/kai/Development/github/NGSovleCoding/data/tw_tuningfork",
    use_input_for_wn=True,
    manifold_surface=True,
    correct_surface_orientation=True)

