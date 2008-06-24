%define _prefix /opt/lscsoft/glue
%define _sysconfdir %{_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_prefix}')")}


Name: 		glue
Summary:	The Grid LSC User Environment
Version:	1.18
Release:	1.lscsoft
License:	None
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	python
BuildRequires:  python-devel
Prefix:         %{_prefix}

%description
Glue (Grid LSC User Environment) is a suite of python modules and programs to
allow users to run LSC codes on the grid.

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

echo  %{name} %{version} [06/19/08]  > %{_prefix}/glue-version

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{python_sitearch}/glue/
%{_prefix}/bin/
%{_prefix}/etc/

%changelog
* Tue Jun 24 2008 Ping Wei <piwei@syr.edu>
- Build for glue 1.18-1

* Fri Jun 20 2008 Xavier Amador <xavier.amador@gravity.phys.uwm.edu>
- added 1 line of code to create "glue-version" file with version info
- under request of Filippo Grimaldi, MIT

* Wed Jun 19 2008 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.17

* Fri Nov 04 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Build for glue 1.6

* Tue Aug 23 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Initial build for glue 1.0
