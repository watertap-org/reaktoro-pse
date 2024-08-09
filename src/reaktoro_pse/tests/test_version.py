import reaktoro_pse


def test_valid_version():
    assert hasattr(reaktoro_pse, "__version__")
    assert isinstance(reaktoro_pse.__version__, str)
