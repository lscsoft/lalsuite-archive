%define _prefix /opt/lscsoft/glue
%define _sysconfdir %{_prefix}/etc
%define _docdir %{_datadir}/doc
%{!?python_sitearch: %define python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1,prefix='%{_prefix}')")}


Name: 		glue
Summary:	The Grid LSC User Environment
Version:	1.25
Release:	1.lscsoft
License:	None
Group:		Development/Libraries
Source:		%{name}-%{version}.tar.gz
Url:		http://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html
BuildRoot:	%{_tmppath}/%{name}-%{version}-root
Requires:	python python-cjson m2crypto
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


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{python_sitearch}/glue/
%{_prefix}/bin/
%{_prefix}/etc/

%changelog
* Wed Jul 01 2009 Duncan Brown <dabrown@physics.syr.edu>
- Pre S6 release of glue

* Wed Jun 24 2009 Duncan Brown <dabrown@physics.syr.edu>
- Post E14 release of glue

* Tue Jun 11 2009 Duncan Brown <dabrown@physics.syr.edu>
- Allow segment tools to see multiple ifos

* Tue Jun 10 2009 Duncan Brown <dabrown@physics.syr.edu>
- Restored LSCdataFindcheck and fixed debian control files

* Tue Jun 09 2009 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.19-1

* Tue Jun 24 2008 Ping Wei <piwei@syr.edu>
- Build for glue 1.18-1

* Wed Jun 19 2008 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.17

* Fri Nov 04 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Build for glue 1.6

* Tue Aug 23 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Initial build for glue 1.0
