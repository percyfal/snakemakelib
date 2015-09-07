# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
index = lambda ref, *args, **kw: ref

def test_index():
    assert "foo.fa" == index("foo.fa", "bar")

def test_cloudbiolinux_import():
    from snakemakelib.bio.ngs.db import index
    assert "foo.fa" == index("foo.fa", "bar")
    try:
        import snakemakelib.bio.ngs.db.cloudbiolinux
    except:
        print("couldn't import snakemakelib.bio.ngs.db.cloudbiolinux")
    from snakemakelib.bio.ngs.db import index
    assert "../bar/foo.fa" == snakemakelib.bio.ngs.db.index("foo.fa", "bar")
    assert "../bar/foo.fa" == index("foo.fa", "bar")
    assert "../bar/foo.fa" == snakemakelib.bio.ngs.db.cloudbiolinux.index("foo.fa", "bar")
