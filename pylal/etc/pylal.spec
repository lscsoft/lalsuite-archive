%define _pylal_prefix /usr
%define _sysconfdir %{_pylal_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?pylal_python_sitearch: %define pylal_python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_pylal_prefix}')")}
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1)")}

Name: 		pylal
Summary:	Python LIGO Algorithm Library
Version:	0.1
Release:	1.lscsoft
License:	See file LICENSE
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/pylal.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	python glue glue-common glue-segments lal lalmetaio lalframe lalburst lalsimulation numpy scipy python-matplotlib
BuildRequires:  python-devel lal-devel lalmetaio-devel lalframe-devel lalburst-devel lalsimulation-devel numpy
Prefix:         %{_pylal_prefix}
%description
The PyLAL package is best described as the Python LIGO Algorithm Library. It was originally a Python wrapping of parts of the LAL library, and although it continues to provide that function it has aquired a large collection of independent code of its own so that it is no longer exclusively a Python interface to LAL.
In this package you will find convenience code to assist with manipulating XML documents using the glue.ligolw I/O library, you will find a wrapping to libframe to enable GWF frame-file reading, you will find binning and smoothing code, and you will find (partial) wrappings of LAL's burstsearch, date, inject, tools, and window packages. Additionally, you will find most, if not all, of the inspiral pipeline's follow-up and summary tools, and several burst-related trigger post-production tools.

%prep
%setup

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install -O1 \
        --skip-build \
        --root=%{buildroot} \
        --prefix=%{_pylal_prefix} \
        --record=INSTALLED_FILES

%clean
rm -rf $RPM_BUILD_ROOT

%files -f INSTALLED_FILES
%defattr(-,root,root)
