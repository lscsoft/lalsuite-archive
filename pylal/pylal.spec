%define _prefix /opt/lscsoft/pylal
%define _sysconfdir %{_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_prefix}')")}


Name: 		pylal
Summary:	LSC Algorithm Library for Python
Version:	0.1
Release:	1.lscsoft
License:	GPL
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/pylal.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	glue lal lalburst libframe
BuildRequires:  lal-devel lalburst-devel libframe-devel python-devel
Prefix:         %{_prefix}

%description
pylal is a collection of python programs and modules for gravitational-wave
data analysis.  pylal modules are a mixture of pure python code and python
bindings of LAL and related libraries.

%prep
%setup 

%build
CFLAGS="%{optflags}" %{__python} setup.py build

%install
rm -rf %{buildroot}
%{__python} setup.py install -O1 \
        --skip-build \
        --root=%{buildroot} \
        --prefix=%{_prefix}


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{python_sitearch}/pylal/
%{_prefix}/bin/
%{_prefix}/etc/
%{_prefix}/var/

%changelog
* Wed Nov 30 2009 Kipp Cannon <kipp.cannon@ligo.org>
- First release of pylal package
