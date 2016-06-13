%define name ligo-lars
%define version 1.1.0
%define unmangled_version 1.1.0
%define release 1

Summary: LIGO Archiving Service
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{unmangled_version}.tar.gz
License: GPL
Group: Development/Libraries
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Prefix: %{_prefix}
BuildArch: noarch
Vendor: Brian Moe <brian.moe@ligo.org>
Requires: ligo-common glue-segments m2crypto
Url: http://www.lsc-group.phys.uwm.edu/daswg/lars.html

%description
Lars is a protype for an archival service for LIGO

%prep
%setup -n %{name}-%{unmangled_version}

%build
python setup.py build

%install
python setup.py install --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES

%clean
rm -rf $RPM_BUILD_ROOT

%files -f INSTALLED_FILES
%defattr(-,root,root)
%exclude %{python_sitelib}/ligo/lars/*pyo
