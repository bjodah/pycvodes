def get_libs(config=None):
    if config is None:
        from . import config
    return (
        "sundials_nvecserial,sundials_cvodes,sundials_sunlinsolspgmr,sundials_sunlinsolspbcgs,"
        "sundials_sunlinsolsptfqmr,sundials_sunmatrixdense,sundials_sunmatrixband" + (
            ',sundials_sunlinsollapackdense,sundials_sunlinsollapackband' if config["LAPACK"]
            else ',sundials_sunlinsoldense,sundials_sunlinsolband'
        ) + (
            ",sundials_sunlinsolklu" if config["KLU"]
            else ""
        )
    )


def get_libs_linkline(config=None):
    libs = get_libs(config)
    if libs:
        return " ".join(["-l%s" % lib for lib in libs.split(",")])
    else:
        return ""


def print_libs_linkline(config=None):
    print(get_libs_linkline(config))
