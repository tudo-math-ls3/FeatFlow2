!##############################################################################
!# ****************************************************************************
!# <name> preproctest </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This simple program contains several typical tests for the
!# preprocessor. Preprocessing is supported natively by some Fortran
!# compilers (using fpp) but others do not provide a dedicated Fortran
!# preprocessor so that the C preprocessor cpp has to be used.  Since
!# cpp does not know about special Fortran syntax, e.g., the maximum
!# length of 132 characters per line and the line continuation symbol
!# '&', cpp is embedded into the f90cpp script which pre-preprocesses
!# the source code passed to cpp and post-preprocesses the outcome in
!# order to handle Fortran specific constructions. Since this is not a
!# trivial task, this test suite contains most of the typical
!# constructions that have to be handled correctly by f90cpp.
!#
!# Moreover, the Featflow2 source code should compile natively under
!# Windows which means that no configure-based build system is
!# available, and thus, the source code must comply with the standard
!# Fortran preprocessor fpp.
!#
!# </purpose>
!##############################################################################

program preproctest

#include "preproctest.h"

  implicit none

  real, dimension(8) :: AVeryLongVariableName = 1
  real, dimension(8) :: YetAnotherVeryLongVariableName = 1

  ! Test #1: Call to macro with very long variable names so that the maximum
  !          length of 132 characters per line is exceeded. There exist the
  !          following remedies to this violation of the Fortran standard:
  !
  !          (a) lines are truncated manually and continued in the next 
  !              line using the keyword 'mynewline' (written in upper case).
  !              This keyword is recognized by f90cpp and converted into '&\n'
  !              after the file has been preprocessed by the preprocessor.
  !              This is the default way if the script f90cpp is used.
  !
  !          (b) the Fortran preprocessor truncates long lines and inserts
  !              line continuation symbols '&' automatically.
  !              This is the case if the SUN Fortran preprocessor is used.
  !
  !          (c) the Fortran compiler does support lines with mode then 132
  !              characters, and therefore, the Fortran preprocessor does
  !              not take special care about extra long lines.
  !              This is the case if the Intel Fortran preprocessor is used.

  write(*,*) "Test #1: Call to macro with very long variable names"
  write(*,*) euclidnorm8(AVeryLongVariableName), sqrt(sum(AVeryLongVariableName**2))
  write(*,*)

  ! Test #2: Call to a macro with arguments over multiple lines.
  !          The maximum length of the macro call would be long and therefore
  !          difficult to read in the source code so that it is truncated and
  !          continued in the next line. Note that no 'mynewline' is required
  !          since we do not want to enforce a linebreak in the preprocessed
  !          code but just extend a macro call over multiple lines.
  !          Depending on the preprocessor the line is either merged into a
  !          singe line before macro expansion or kept as is.

  write(*,*) "Test #2: Call to a macro over multiple lines"
  write(*,*) combine(\
    euclidnorm8(AVeryLongVariableName),\
    euclidnorm8(YetAnotherVeryLongVariableName))
  write(*,*)

  ! Test #3: Call to a macro inside a list of other arguments.
  !          This test extends Test #2 in the following sense. The preprocessor
  !          must not add extra lines after the expanded macro call because then
  !          the Fortran compiler may complain about unbalance parantheses.
  !          For this reasons, newer versions of the GNU cpp must be envoked
  !          the '--traditional' flag. Otherwise, the following test fails.

  write(*,*) "Test #3: Call to a macro inside other arguments"
  call writeTwoArguments(combine(\
  euclidnorm8(AVeryLongVariableName),\
  euclidnorm8(YetAnotherVeryLongVariableName)),&
      sqrt(sum(AVeryLongVariableName**2))+&
      sqrt(sum(YetAnotherVeryLongVariableName**2)))
  write(*,*)


  ! Test #4: The preprocessor should not manipulate comments which include 
  !          characters such as backslashes which would normally be interpreted
  !          as line continuation symbols by the C preprocessor. This special
  !          case is handled by perl command (4) in script f90cpp.
  !
  !          It applies to ASCII arts of the form
  !
  !          / a_11  a_12  a_13  a_14 \
  !          | a_21  a_22  a_23  a_24 |
  !          | a_31  a_32  a_33  a_34 |
  !          \ a_41  a_42  a_43  a_44 /
  !
  !          as well as to LaTeX forced line breaks in comments, e.g.
  !
  !          $$ M = left(\begin{array}{cccc}
  !                      a_{11} & a_{12} & a_{13} & a_{14}\\
  !                      a_{21} & a_{22} & a_{23} & a_{24}\\
  !                      a_{31} & a_{32} & a_{33} & a_{34}\\
  !                      a_{41} & a_{42} & a_{43} & a_{44}
  !                 \end{array}\right) $$
  
  ! Remark: In the write statements below it is necessary to have a white space
  !         between the backslash and the quotation mark. This has nothing to
  !         to with the preprocessor but could be misinterpreted by some Fortran
  !         compilers (e.g., PGI Fortran 90 Compiler version 9.0-4). Not yet sure
  !         if this is a bug in the PGI compiler or not valid according to the
  !         Fortran 90 standard.

  write(*,*) "Test #4: Correct treatment of special characters in comments"
  write(*,*) "/ a_11  a_12  a_13  a_14 \ "
  write(*,*) "| a_21  a_22  a_23  a_24 |"
  write(*,*) "| a_31  a_32  a_33  a_34 |"
  write(*,*) "\ a_41  a_42  a_43  a_44 /"
  write(*,*)


  ! Test #5: The preprocessor must not interpred Fortran 90 string 
  !          concatination operators ('//'). This special case is
  !          handled by perl command (2) in script f90cpp.

  write(*,*) "Test #5: Correct treatment of Fortran 90 string concatenation"
  write(*,*) "Foo"//"bar"
  write(*,*)

  
  ! Test #6: Stringification of a macro argument. This feature is supported
  !          by GNU cpp if NOT envoked in traditional mode as well as Fortran
  !          preprocessors (tested with Intel and SUN fpp). However, newer
  !          versions of the GNU cpp preprocessor must be envoked with the
  !          '--traditional flag' enabled in order to handle Test #3 correctly.
  !          Therefore, stringification of macro arguments is not a feature
  !          that has to be supported on all platforms, and thus, it should
  !          not be used in the Featflow2 source code.

!!$  write(*,*) "Test #6: Stringification of macro arguments"
!!$  write(*,*) stringify(print)
!!$  write(*,*)

  
  ! Test #7: Concatenation of two arguments. This feature is supported by 
  !          GNU cpp if NOT envoked in traditional mode as well as the
  !          Intel Fortran preprocessor. It is not supported by the SUN
  !          Fortran preprocessor and the GNU cpp envoked in traditional
  !          mode. Therefore, concatenation of strings is not a feature
  !          that has to be supported on all platforms, and thus, it should
  !          not be used in the Featflow2 source code.

!!$  write(*,*) "Test #7: Concatenation of two macro arguments"
!!$  write(*,*) concatenate(foo,bar)
!!$  write(*,*)

contains

  subroutine writeTwoArguments(a,b)
    real, intent(in) :: a,b

    write(*,*) a,b

  end subroutine writeTwoArguments

end program preproctest

