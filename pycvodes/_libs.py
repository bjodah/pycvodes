def _sundials_major_version():
    import os
    import re
    cxxflags = os.environ.get('CXXFLAGS', '')
    sundials_inc = os.environ.get('SUNDIALS_INC', '')
    inc_dirs = [sundials_inc] if sundials_inc else []
    # Parse -isystem and -I flags, handling both "-isystem/path" and "-isystem /path"
    parts = cxxflags.split()
    i = 0
    while i < len(parts):
        part = parts[i]
        if part in ('-isystem', '-I') and i + 1 < len(parts):
            inc_dirs.append(parts[i + 1])
            i += 2
        elif part.startswith('-isystem'):
            inc_dirs.append(part[len('-isystem'):])
            i += 1
        elif part.startswith('-I') and len(part) > 2:
            inc_dirs.append(part[2:])
            i += 1
        else:
            i += 1
    for inc_dir in inc_dirs:
        cfg = os.path.join(inc_dir, 'sundials', 'sundials_config.h')
        if os.path.isfile(cfg):
            txt = open(cfg).read()
            m = re.search(r'#define\s+SUNDIALS_VERSION_MAJOR\s+(\d+)', txt)
            if m:
                return int(m.group(1))
    return None


def get_libs(config=None):
    if config is None:
        from . import config
    sundials_major = _sundials_major_version()
    core_lib = ",sundials_core" if (sundials_major is not None and sundials_major >= 7) else ""
    return (
        "sundials_nvecserial,sundials_cvodes,sundials_sunlinsolspgmr,sundials_sunlinsolspbcgs,"
        "sundials_sunlinsolsptfqmr,sundials_sunmatrixdense,sundials_sunmatrixband,sundials_sunmatrixsparse" +
        core_lib + (
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
