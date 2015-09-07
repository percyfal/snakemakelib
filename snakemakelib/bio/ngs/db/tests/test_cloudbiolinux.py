# Copyright (C) 2015 by Per Unneberg
# pylint: disable=R0904
import pytest
index = lambda ref, *args, **kw: ref

def test_index():
    assert "foo.fa" == index("foo.fa", "bar")

def test_cloudbiolinux_import():
    """Test that importing cloudbiolinux changes db.index"""
    from snakemakelib.bio.ngs.db import index
    assert "foo.fa" == index("foo.fa", "bar", index="foo.fa")
    import snakemakelib.bio.ngs.db.cloudbiolinux
    from snakemakelib.bio.ngs.db import index
    assert "../bar/foo.fa" == snakemakelib.bio.ngs.db.index("foo.fa", "bar")
    assert "../bar/foo.fa" == index("foo.fa", "bar")
    assert "../bar/foo.fa" == snakemakelib.bio.ngs.db.cloudbiolinux.index("foo.fa", "bar")

def test_ref_wo_build():
    """Test reference function return value when build is missing"""
    from snakemakelib.bio.ngs.db.cloudbiolinux import ref
    cfg = {'bio.ngs.settings' : {'db' : {'build' : None, 'ref':"foo/bar"}}}
    assert ref("foobar", cfg['bio.ngs.settings']['db']) == "foo/foobar"

def test_ref_wo_build_config():
    """Test reference function return value when build present but build_config is not"""
    from snakemakelib.bio.ngs.db.cloudbiolinux import ref
    cfg = {'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : "foo/bar"}}}
    assert ref("foobar", cfg['bio.ngs.settings']['db']) ==  "foo/foobar"


def test_index_missing_application():
    """Test calling index without applying application"""
    from snakemakelib.bio.ngs.db.cloudbiolinux import index
    config = {'bio.ngs.settings' : {'db' : {'build' : 'hg19', 'build_config' : None, 'ref' : "foo.fa"}}, 'bio.ngs.align.bwa': {'index' : index}}
    bwa_cfg = config['bio.ngs.align.bwa']
    with pytest.raises(TypeError):
        index(ref=config['bio.ngs.settings']['db']['ref'], index=bwa_cfg['index'])
