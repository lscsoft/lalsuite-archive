# SWIG configuration
# Author: Karl Wette, 2011
#
# serial 6

# basic version string comparison
# can only handle numeric versions separated by periods
AC_DEFUN([LALSUITE_VERSION_COMPARE],[
  vcmp_awkprog='{n = 0; while(++n <= NF) { if($n > 99) printf "ERROR"; printf "%02u", $n; } }'
  vcmp_v1=[`echo $1 | ${SED} '/[^0-9.]/d' | ${AWK} -F . "${vcmp_awkprog}" | ${SED} '/ERROR/d'`]
  vcmp_v2=[`echo $2 | ${SED} '/[^0-9.]/d' | ${AWK} -F . "${vcmp_awkprog}" | ${SED} '/ERROR/d'`]
  AS_IF([test x${vcmp_v1} = x],[AC_MSG_ERROR([could not parse version string '$1'])])
  AS_IF([test x${vcmp_v2} = x],[AC_MSG_ERROR([could not parse version string '$2'])])
  AS_IF([test ${vcmp_v1} -lt ${vcmp_v2}],[$3])
  AS_IF([test ${vcmp_v1} -eq ${vcmp_v2}],[$4])
  AS_IF([test ${vcmp_v1} -gt ${vcmp_v2}],[$5])
])

# SWIG setup and configuration
AC_DEFUN([LALSUITE_WITH_SWIG],[

  # minimum required SWIG version
  SWIG_MIN_VERSION=1.3.40

  # save and clear CPPFLAGS and LIBS
  swig_CPPFLAGS=${CPPFLAGS}
  swig_LIBS=${LIBS}
  CPPFLAGS=
  LIBS=

  # check for sed
  AC_PROG_SED

  # check for awk
  AC_PROG_AWK

  # check for MKDIR_P
  m4_ifdef([AC_PROG_MKDIR_P],[],[
    MKDIR_P='$(INSTALL) -d'
    AC_SUBST(MKDIR_P)
  ])

  # ./configure command line option
  AC_ARG_WITH(
    [swig],
    AC_HELP_STRING(
      [--with-swig=LANG...],
      [generate SWIG wrappings for scripting languages LANG]
    ),
    [swig_langs="${with_swig}"],   # --with-swig=...
    [swig_langs=no]                # --without-swig / no option
  )

  # replace any delimiting commas with spaces
  swig_langs=[`echo ${swig_langs} | ${SED} 's|,| |g'`]

  # are we are binding LAL itself, or one of the other LAL libraries?
  swig_islal="test x${PACKAGE_NAME} = xlal"

  # are we (not) in debugging mode?
  swiglal_ndebug="test x${enable_debug} = xno"

  # common SWIG interface headers (with LAL only)
  AS_IF([${swig_islal}],[
    SWIG_HEADERS=
    AC_SUBST(SWIG_HEADERS)
  ])

  # list of SWIG wrappings which will be built
  SWIG_WRAPPINGS=

  # iterate over SWIG target scripting languages
  swig_any_lang_found=
  for swig_lang in ${swig_langs}; do
    swig_lang_found=

    # configure this language, or fail if language was not found
    LALSUITE_SWIG_LANGUAGES
    AS_IF([test "x${swig_lang_found}" = xyes],[
      swig_any_lang_found=yes
    ],[
      AS_IF([test "x${swig_lang}" != xno],[
        AC_MSG_ERROR([no SWIG support for LANG='${swig_lang}'])
      ])
    ])

  done

  # check if any language was found
  AS_IF([test x${swig_any_lang_found} = xyes],[
    swig_build_cond="test 1 == 1"
  ],[
    swig_build_cond="test 1 == 0"
  ])
  AM_CONDITIONAL(SWIG_BUILD,[${swig_build_cond}])
  AS_IF([${swig_build_cond}],[

    # common SWIG language build makefile
    SWIG_COMMON_MK="${srcdir}/../lal/swig/swig-common.mk"
    AC_SUBST_FILE(SWIG_COMMON_MK)

    # SWIG makefile for generating header lists
    SWIG_HEADER_MK="${srcdir}/../lal/swig/swig-header.mk"
    AC_SUBST_FILE(SWIG_HEADER_MK)

    # symbols to define when generating SWIG wrapping code
    SWIG_SWIG_DEFINES=
    AC_SUBST(SWIG_SWIG_DEFINES)

    # symbols to define when compiling SWIG wrapping code
    SWIG_CXX_DEFINES=
    AC_SUBST(SWIG_CXX_DEFINES)

    # are we are binding LAL itself, or one of the other LAL libraries?
    AS_IF([${swig_islal}],[
      SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_IS_LAL"
    ])

    # are we (not) in debugging mode?
    AS_IF([${swiglal_ndebug}],[
      SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_NDEBUG"
      SWIG_CXX_DEFINES="${SWIG_CXX_DEFINES} SWIGLAL_NDEBUG"
    ])

    # try to figure out the underlying type of int64_t
    AC_CHECK_HEADERS([stdint.h],[],[
      AC_MSG_ERROR([could not find 'stdint.h'])
    ])
    AC_MSG_CHECKING([underlying type of int64_t])
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],[
          int64_t i64 = 0; const long int *pli = &i64;
        ])
      ],[
        AC_MSG_RESULT([long int])
        swig_wordsize=SWIGWORDSIZE64
      ],[
        AC_COMPILE_IFELSE(
          [
            AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],[
              int64_t i64 = 0; const long long int *plli = &i64;
            ])
          ],[
            AC_MSG_RESULT([long long int])
            swig_wordsize=
          ],[
            AC_MSG_ERROR([could not determine underlying type of int64_t])
          ]
        )
      ]
    )
    AC_LANG_POP([C++])
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} ${swig_wordsize}"

    # ensure that all LAL library modules share type information
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIG_TYPE_TABLE=swiglaltypetable"

    # make SWIG use C++ casts
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIG_CPLUSPLUS_CAST"

    # define C99 constant and limit macros
    SWIG_CXX_DEFINES="${SWIG_CXX_DEFINES} __STDC_CONSTANT_MACROS __STDC_LIMIT_MACROS"               

    # common SWIG interface headers (with LAL only)
    AS_IF([${swig_islal}],[
      SWIG_HEADERS="${SWIG_HEADERS} swiglal-common.i swiglal-gsl.i"
    ])

    # path SWIG should look in for header files:
    # keep any -I options in CPPFLAGS, without the -I prefix
    SWIG_INCLPATH=[`for n in ${swig_CPPFLAGS}; do echo $n | ${SED} 's|^-I||p;d'; done`]
    SWIG_INCLPATH=[`echo ${SWIG_INCLPATH}`]   # get rid of newlines
    AC_SUBST(SWIG_INCLPATH)

    # path SWIG should look in for libraries:
    # keep any -L options in _LIB variables, without the -L prefix
    # keep any "lib*.la" files, replace filename with $objdir;
    swig_regex=['s|^-L||p;s|lib[^/][^/]*\.la|'"${objdir}"'|p;d']
    SWIG_LIBPATH="${LAL_LIBS} ${LALSUPPORT_LIBS} ${swig_LIBS}"
    SWIG_LIBPATH=[`for n in ${SWIG_LIBPATH}; do echo $n | ${SED} "${swig_regex}"; done`]
    SWIG_LIBPATH=[`echo ${SWIG_LIBPATH}`]   # get rid of newlines
    # add pre-install locations for lal, lalsupport, and lal* libraries
    SWIG_LIBPATH="${SWIG_LIBPATH} \$(top_builddir)/lib/${objdir}"
    SWIG_LIBPATH="${SWIG_LIBPATH} \$(top_builddir)/packages/support/src/${objdir}"
    SWIG_LIBPATH="${SWIG_LIBPATH} \$(top_builddir)/src/${objdir}"
    AC_SUBST(SWIG_LIBPATH)

  ],[

    # if no SWIG languages were found
    SWIG_WRAPPINGS="NONE"
    SWIG_COMMON_MK="/dev/null"
    SWIG_HEADER_MK="/dev/null"

  ])

  # restore CPPFLAGS and LIBS
  CPPFLAGS=${swig_CPPFLAGS}
  LIBS=${swig_LIBS}

])

# tell the SWIG wrappings to use some feature
AC_DEFUN([LALSUITE_SWIG_USE],[
  SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_USE_$1"
])

# SWIG language configuration
AC_DEFUN([LALSUITE_SWIG_LANG],[

  # check whether to configure this language
  AS_IF([test x${swig_lang} = x$1 || test x${swig_lang} = xyes],[
    swig_build_lang_cond="test 1 == 1"
  ],[
    swig_build_lang_cond="test 1 == 0"
  ])
  m4_pushdef([swigcondname],translit([swig_build_$1],[a-z],[A-Z]))
  AM_CONDITIONAL(swigcondname,[${swig_build_lang_cond}])
  AS_IF([${swig_build_lang_cond}],[

    # first-time general configuration
    AS_IF([test "x${swig_any_lang_found}" = x],[

      # check for swig binary
      AC_PATH_PROGS(SWIG,[swig],[])
      AS_IF([test "x${SWIG}" = x],[
        AC_MSG_ERROR([could not find 'swig' in path])
      ])

      # check for swig version
      AC_MSG_CHECKING([for swig version])
      swig_regex=['s|^.*SWIG\s\s*[Vv]ersion\s\s*\([0-9.][0-9.]*\).*$|\1|p;d']
      swig_version=[`${SWIG} -version | ${SED} "${swig_regex}"`]
      AS_IF([test "x${swig_version}" = x],[
        AC_MSG_ERROR([could not determine swig version])
      ])
      AC_MSG_RESULT([${swig_version}])

      # check if swig version is newer than required
      LALSUITE_VERSION_COMPARE([${SWIG_MIN_VERSION}],[${swig_version}],[],[],[
        AC_MSG_ERROR([require swig version >= ${SWIG_MIN_VERSION}])
      ])

    ])

    # language was found, so add to list
    swig_lang_found=yes
    AS_IF([test "x${SWIG_WRAPPINGS}" = x],[
      SWIG_WRAPPINGS="$1"
    ],[
      SWIG_WRAPPINGS="${SWIG_WRAPPINGS}, $1"
    ])

    # configure language
    $2

    # language-specific SWIG interface header (with LAL only)
    AS_IF([${swig_islal}],[
      SWIG_HEADERS="${SWIG_HEADERS} $1/swiglal-$1.i"
    ])

  ])
  m4_popdef([swigcondname])

])

# SWIG languages
AC_DEFUN([LALSUITE_SWIG_LANGUAGES],[
  LALSUITE_SWIG_LANG_OCTAVE
  LALSUITE_SWIG_LANG_PYTHON
])

# SWIG octave configuration
AC_DEFUN([LALSUITE_SWIG_LANG_OCTAVE],[
  LALSUITE_SWIG_LANG([octave],[

    # minimum required octave version
    OCTAVE_MIN_VERSION=3.2.4

    # check for octave-config binary
    AC_PATH_PROGS(OCTAVE_CONFIG,[octave-config],[])
    AS_IF([test "x${OCTAVE_CONFIG}" = x],[
      AC_MSG_ERROR([could not find 'octave-config' in path])
    ])

    # check for corresponding octave binary
    AC_MSG_CHECKING([for octave])
    OCTAVE=`${OCTAVE_CONFIG} -p BINDIR`/octave
    AS_IF([test -f "${OCTAVE}" && test -x "${OCTAVE}"],[],[
      AC_MSG_ERROR([could not find 'octave' in path])
    ])
    AC_MSG_RESULT([${OCTAVE}])
    AC_SUBST(OCTAVE)
    # add flags for silence and environment-independence
    OCTAVE="${OCTAVE} -qfH"

    # check for octave version
    AC_MSG_CHECKING([for octave version])
    octave_version=`${OCTAVE_CONFIG} --version`
    AS_IF([test "x${octave_version}" = x],[
      AC_MSG_ERROR([could not determine octave version])
    ])
    AC_MSG_RESULT([${octave_version}])

    # check if octave version is newer than required
    LALSUITE_VERSION_COMPARE([${OCTAVE_MIN_VERSION}],[${octave_version}],[],[],[
      AC_MSG_ERROR([require octave version >= ${OCTAVE_MIN_VERSION}])
    ])

    # determine where to install .oct files:
    # take site .oct install dir given by octave-config,
    # and strip off prefix; thus, if LAL is installed in
    # the same directory as octave, .oct files will be
    # found by octave without having to add to OCTAVE_PATH
    octave_prefix=[`${OCTAVE_CONFIG} -p PREFIX | ${SED} 's|/*$||'`]
    AC_MSG_CHECKING([for octave .oct installation directory])
    octave_localoctfiledir=[`${OCTAVE_CONFIG} -p LOCALOCTFILEDIR | ${SED} 's|/*$||'`]
    OCTAVE_OCTFILEDIR=[`echo ${octave_localoctfiledir} | ${SED} "s|^${octave_prefix}/||"`]
    AS_IF([test -n "`echo ${OCTAVE_OCTFILEDIR} | ${SED} -n '\|^/|p'`"],[
      AC_MSG_ERROR([could not build relative path from '${OCTAVE_OCTFILEDIR}'])
    ])
    AC_MSG_RESULT([${OCTAVE_OCTFILEDIR}])
    AC_SUBST(OCTAVE_OCTFILEDIR)

  ])
])

# SWIG python configuration
AC_DEFUN([LALSUITE_SWIG_LANG_PYTHON],[
  LALSUITE_SWIG_LANG([python],[

    # check for python
    AM_PATH_PYTHON([2.4])

    # check for numpy
    AC_MSG_CHECKING([for numpy])
    ${PYTHON} -c "import numpy" 2>/dev/null
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not import numpy])
    ])
    AC_MSG_RESULT([yes])

  ])
])
