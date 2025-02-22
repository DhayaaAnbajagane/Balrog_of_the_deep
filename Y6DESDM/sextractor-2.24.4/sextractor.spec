#
#				sextractor.spec.in
#
# Process this file with autoconf to generate an RPM .spec packaging script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SExtractor
#
#	Copyright:		(C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SExtractor is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SExtractor is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		19/06/2017
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define name sextractor
%define version 2.24.4
%define release 1%{?dist}
%undefine _missing_build_ids_terminate_build

Summary: Extract catalogs of sources from astronomical images
Name: %{name}
Version: %{version}
Release: %{release}
Source0: ftp://ftp.iap.fr/pub/from_users/bertin/sextractor/%{name}-%{version}.tar.gz
URL: http://astromatic.net/software/%{name}
License: GPL v3+
Group: Sciences/Astronomy
BuildRoot: %{_tmppath}/%{name}-buildroot
BuildRequires: pkgconfig
BuildRequires: fftw-devel >= 3.1
BuildRequires: atlas-devel >= 3.6.0

%description
Extract catalogs of sources from astronomical images

%prep
%setup -q

%build
if test "$USE_BEST"; then
%configure --enable-mkl --enable-auto-flags --enable-best-link
elif test "$USE_ICC"; then
%configure --enable-icc
else
%configure
fi
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc config AUTHORS BUGS ChangeLog COPYRIGHT HISTORY INSTALL LICENSE README.md THANKS
%{_bindir}/sex
%{_bindir}/ldactoasc
%{_mandir}/man1/sex.1*
%{_mandir}/manx/sex.x*
%{_datadir}/sextractor

%changelog
* Wed Apr 27 2022 AstrOmatic <astromatic@astromatic.net>
- Automatic RPM rebuild
* Tue May 13 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3
* Fri Apr 04 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3b4
* Wed Mar 05 2003 Emmanuel Bertin <bertin@iap.fr>
- RPM build for V2.3b3
* Fri Feb 07 2003 Emmanuel Bertin <bertin@iap.fr>
- Second RPM build
* Fri Jan 24 2003 Emmanuel Bertin <bertin@iap.fr>
- Second RPM build
* Sun Dec 15 2002 Emmanuel Bertin <bertin@iap.fr>
- First RPM build

# end of file
