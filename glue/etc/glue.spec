%define _prefix /opt/lscsoft/glue
%define _sysconfdir %{_prefix}/etc
%define _docdir %{_datadir}/doc

Name: 		glue
Summary:	The Grid LSC User Environment
Version:	1.17
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
python setup.py build

%install
rm -rf $RPM_BUILD_ROOT
python setup.py install --root=${RPM_BUILD_ROOT} --prefix=/opt/lscsoft/glue --record=INSTALLED_FILES

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{prefix}


%changelog
* Wed Jun 11 2008 Duncan Brown <dabrown@physics.syr.edu>
- Build for glue 1.17

* Fri Nov 04 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Build for glue 1.6

* Tue Aug 23 2005 Duncan Brown <dbrown@ligo.caltech.edu>
- Initial build for glue 1.0
