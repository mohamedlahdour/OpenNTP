! Takes Fortran 77 code in standard format and makes some changes to produce
! free-format Fortran 90 code.
! N.B. It expects STANDARD F77 code.   Non-standard extensions such as
! DO .. END DO (i.e. no label) or in-line comments may cause havoc!

! Changes included are:
!     C or c in column 1 replaced with !
!     Continuation denoted by a character in column 6 replaced with & at the
!         end of the previous line.
!     Indenting of code for DO-loops and IF blocks.
!     END of program unit replaced by END SUBROUTINE (/PROGRAM/FUNCTION) name
!     Fortran `keywords' are in upper case, all other words other than those
!         in character strings are converted to lower case.
!     .LT., .EQ., etc. replaced with <, ==, etc.
!     Labels removed from DO loops; all of which will end with END DO.
!         If labels are not referenced, they are removed.
!     Short continued lines are adjoined to the previous line.
!     ENDIF, ELSEIF & GOTO split into separate words.
!     3-way arithmetic IF constructs are converted to IF .. ELSE IF form.
!     Embedded blanks are removed from numbers in DATA statements.
!     INTENT declarations are added for dummy arguments.
!     Some GO TOs are converted to CYCLE or EXIT.
!     Converts CHARACTER * to CHARACTER (LEN=xx) ::.
!     Converts computed GO TOs to SELECT CASE.

! To be done:
!     DATA statements to be replaced by assignments on the declaration line.
!     IMPLICIT NONE statements to be included.
!     Declaration of types of unlisted variables.
!     Functions to be converted to ELF90 form, i.e. REAL FUNCTION XYZ(arg)
!         converted to FUNCTION xyz(arg) RESULT(fn_val).

! Known problems
!     Cannot handle character strings or names broken at the end of lines.
!     No attempt to convert BLOCKDATA, COMMON or EQUIVALENCE.
!     Does not convert Hollerith strings, e.g. 31HTHIS IS A COMMENT ...
!     May do the wrong thing if variable names start with IF or end with DO.
!     INTENTs are sometimes wrong.  In particular, INTENT(IN) arguments are
!         often shown as INTENT(IN OUT).
!     Cannot handle comment lines in the middle of continued instructions.
!     Can handle 'character*(*) str' but not 'character str*(*)'.

! The default extension for the name of the input file is `for'; this can be
! over-ruled by giving the full name (e.g. myprog.f77).   The output file name
! will be the input name (and directory) with extension `.f90'.

! Added conversion of `enddo' to END DO - 13 March 1997
! Corrected bug which occurred when an arithmetic IF within a DO-loop involved
!     a jump to the end of the DO-loop - 17 August 1997.

! ELSEIF, ENDIF & ELSEIF were being split into 2 separate words, and then the
!     last letter converted back to lower case - corrected 17 August 1997.
! Corrected bug which occurred when .LT. (or other comparison) had a blank
!     before and/or after, followed on the same line by a text string, followed
!     by a Fortran word such as THEN or GO TO - 8 December 1997.
! Added (LEN=1) after CHARACTER if length not specified - 9 December 1997.
! Embedded blanks are removed from numerical constants in DATA statements.
!     Added 9 December 1997.
! Added INTENTs and TYPE declarations for dummy arguments - 23 December 1997.
! Corrected problem when DO statement contains a comma immediately after DO,
!     and improved the detection of INTENTs when a dummy argument appears in an
!     IF-expression.  Added extra indentation on continuation lines to improve
!     readability - 13 January 1998
! Corrected a bug which could occur when the last type declaration was matched
!     to a dummy variable and the line deleted - 5 June 1998
! Corrected jumps out of inner nested DO loops, and replaced GO TOs, out of
!     DO loops to the next executable line, with EXIT - 8 June 1998
! Added conversion of CHARACTER * to CHARACTER (LEN=xx) ::
!     including CHARACTER*10 a, d, c*50, d   - 21 June 1998.
! Corrected for case of final command of a DO loop which is not CONTINUE and
!     which flows onto the next line - 29 June 1998.
! Added conversion of computed GO TOs to SELECT CASE form, and
!     fixed problem when a CHARACTER dummy variable had '*' as its last
!     dimension - 26 November 1998.
! Fixed problems when the dimensions of a dummy variable involved another
!     dummy variable, e.g. wk(nrow*(ncols+1)) - 25 December 1998
! Added date & time stamp - 27 January 1999
! Finally fixed the problems with CYCLE & EXIT, I hope! - 2 February 1999
! Fixed a problem when a type declaration was continued and the next line
!     declared the type(s) of dummy arguments - 3 February 1999
! Added conversion of PARAMETER statements from PARAMETER (name1=v1, .. )
!     to TYPE1, PARAMETER :: name1=v1  - 8 February 1999
! Added EQV to the list of FORTRAN `words' - 11 February 1999
! Partially corrected problems with the construct:
!     IF (condition) GO TO (10, 20, ..
!       ..., 99), next
!     i.e. with IF & computed GOTO in the same statement (without a THEN), and
!     continued onto the next line.
!     Also changed a DATA statement to a PARAMETER statement to make the code
!     compatible with ELF90 (Thanks to David Ormerod) - 20 May 1999
! Added test for existence of source file.  Program crashed previously if
!     the file was not found - 3 August 1999
! Corrected SUBROUTINE fix_3way_IF so that it does not interpret IFIX (or
!     similar) as an IF - 23 January 2000.
! At last fixed up strings in quotes which flowed from one line to the next
!     - 24 January 2000
! Fixed an error which sometimes caused GOTOs to be wrongly converted to CYCLE
!     - 21 March 2000

! Latest revision - 21 March 2000
! Author - Alan Miller  (amiller @ bigpond.net.au)
! WWW-page: http://users.bigpond.net.au/amiller/


MODULE implicit
! Module to set and reset implicit variable types for use by to_f90.

IMPLICIT NONE
INTEGER, SAVE :: var_type(26) = (/  &
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1 /)
!                a b c d e f g h i j k l m n o p q r s t u v w x y z
CHARACTER (LEN=24), SAVE :: vt(0:7) = (/ 'NO TYPE                 ', &
                                         'REAL                    ', &
                                         'INTEGER                 ', &
                                         'DOUBLE PRECISION        ', &
                                         'LOGICAL                 ', &
                                         'COMPLEX                 ', &
                                         'CHARACTER               ', &
                                         'OTHER TYPE              ' /)

CONTAINS


SUBROUTINE reset_defaults()

var_type(1:8) = 1            ! REAL (A-H)
var_type(9:14) = 2           ! INTEGER (I-N)
var_type(15:26) = 1          ! REAL (O-Z)

RETURN
END SUBROUTINE reset_defaults



SUBROUTINE set_implicit_types(text)
! Read in implicit statement and interpret.

CHARACTER (LEN=*), INTENT(IN OUT) :: text

! Local variables
INTEGER :: ivt, length, start, i, j, pos, left, right
LOGICAL :: first

i = INDEX(text, 'IMPLICIT')
IF (i > 0) text = text(i+8:)
text = ADJUSTL(text)

DO
  IF (text(1:4) == 'NONE') THEN
    var_type = 0
    RETURN
  ELSE IF (text(1:4) == 'REAL') THEN
    ivt = 1
  ELSE IF (text(1:7) == 'INTEGER') THEN
    ivt = 2
  ELSE IF (text(1:24) == 'DOUBLE PRECISION COMPLEX') THEN
    ivt = 7
    vt(7) = 'DOUBLE PRECISION COMPLEX'
  ELSE IF (text(1:16) == 'DOUBLE PRECISION') THEN
    ivt = 3
  ELSE IF (text(1:7) == 'LOGICAL') THEN
    ivt = 4
  ELSE IF (text(1:7) == 'COMPLEX') THEN
    ivt = 5
  ELSE IF (text(1:9) == 'CHARACTER') THEN
    ivt = 6
  ELSE
    ivt = 7
    i = INDEX(text, ' ')
    vt(7) = text(1:i-1)
  END IF

! Interpret the part in brackets, e.g. (a - h, o - z)

  length = LEN_TRIM(text)
  start = 5
  left = INDEX(text(start:length), '(') + start - 1
  IF (left < start) RETURN
  right = INDEX(text(start:length), ')') + start - 1
  IF (right < left) RETURN
                                       ! Interpret text(left+1:right-1)
  first = .TRUE.
  DO pos = left+1, right
    SELECT CASE (text(pos:pos))
      CASE (' ')
        CYCLE
      CASE ('-')
        first = .FALSE.
      CASE (',', ')')
        IF (first) THEN
          var_type(i) = ivt
        ELSE
          var_type(i:j) = ivt
          first = .TRUE.
        END IF
      CASE DEFAULT
        IF (first) THEN
          i = ICHAR(text(pos:pos)) - ICHAR('a') + 1
          IF (i < 1) THEN
            i = ICHAR(text(pos:pos)) - ICHAR('A') + 1
          END IF
        ELSE
          j = ICHAR(text(pos:pos)) - ICHAR('a') + 1
          IF (j < 1) THEN
            j = ICHAR(text(pos:pos)) - ICHAR('A') + 1
          END IF
        END IF
    END SELECT
  END DO

  start = right + 1
  IF (start >= length) RETURN
  text = text(start:length)
  DO
    IF (text(1:1) == ',' .OR. text(1:1) == ' ') THEN
      text = text(2:)
    ELSE
      EXIT
    END IF
  END DO
END DO

RETURN
END SUBROUTINE set_implicit_types



FUNCTION implicit_type(ch) RESULT(vtype)
! Return the variable type given the first character of its name.
! The first character is expected to be lower case, but just in case ..

CHARACTER (LEN=1), INTENT(IN) :: ch
CHARACTER (LEN=24)            :: vtype

! Local variable
INTEGER  :: i, j

i = ICHAR(ch) - ICHAR('a') + 1
IF (i >= 1 .AND. i <= 26) THEN
  j = var_type(i)
  vtype = vt(j)
ELSE
  i = ICHAR(ch) - ICHAR('A') + 1
  IF (i >= 1 .AND. i <= 26) THEN
    j = var_type(i)
    vtype = vt(j)
  ELSE
    vtype = ' '
  END IF
END IF

RETURN
END FUNCTION implicit_type

END MODULE implicit



PROGRAM to_f90
USE implicit
IMPLICIT NONE

TYPE :: code
  CHARACTER (LEN=140)  :: text
  CHARACTER (LEN=  5)  :: label
  TYPE (code), POINTER :: next
END TYPE code

TYPE :: argument
  CHARACTER (LEN=10)       :: name
  INTEGER                  :: intention    ! IN = 1, OUT = 2, IN OUT = 3
  CHARACTER (LEN=24)       :: var_type     ! Room for DOUBLE PRECISION COMPLEX
  INTEGER                  :: dim          ! DIM = 0 for scalars
  CHARACTER (LEN=24)       :: dimensions   ! Not used if DIM = 0
  TYPE (argument), POINTER :: next
END TYPE argument

CHARACTER (LEN=60)       :: f77_name, f90_name
CHARACTER (LEN= 1)       :: tab = CHAR(9), ch
CHARACTER (LEN=50)       :: prog_unit_name = ' ', blank = ' ', case_expr
CHARACTER (LEN= 9)       :: delimiters = ' =+-*/,()'
CHARACTER (LEN=10)       :: numbers = '1234567890'
CHARACTER (LEN= 5)       :: lab
CHARACTER (LEN=30)       :: text, vtype
CHARACTER (LEN=140)      :: statement
CHARACTER (LEN= 8)       :: date
CHARACTER (LEN=10)       :: time
INTEGER                  :: iostatus, pos, count, last, n_marks, pos1(20),  &
                            pos2(20), lab_length, indent, i, i1, i2, length, &
                            numb_arg, i3, i4
TYPE (code), POINTER     :: head, current, tail, last_line, next_line,  &
                            first_decl, last_decl, start_prog_unit,  &
                            end_prog_unit
LOGICAL                  :: asterisk, OK, data_stmnt, first_arg, continuation
TYPE (argument), POINTER :: arg_start, arg, last_arg

INTERFACE
  SUBROUTINE mark_text(text, n_marks, pos1, pos2, continuation)
    IMPLICIT NONE
    CHARACTER (LEN = *), INTENT(IN)  :: text
    INTEGER, INTENT(OUT)             :: n_marks, pos1(:), pos2(:)
    LOGICAL, INTENT(IN)              :: continuation
  END SUBROUTINE mark_text

  SUBROUTINE convert_text(text, n_marks, pos1, pos2)
    IMPLICIT NONE
    CHARACTER (LEN = *), INTENT(IN OUT) :: text
    INTEGER, INTENT(IN)                 :: n_marks
    INTEGER, INTENT(IN OUT)             :: pos1(:), pos2(:)
  END SUBROUTINE convert_text

  SUBROUTINE remove_data_blanks(text)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN OUT) :: text
  END SUBROUTINE remove_data_blanks

  FUNCTION last_char( text ) RESULT(ch)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: text
    CHARACTER (LEN=1)             :: ch
  END FUNCTION last_char

  FUNCTION find_delimited_name (text, name) RESULT(pos)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: text, name
    INTEGER                       :: pos
  END FUNCTION find_delimited_name
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO')' Enter name of Fortran source file: '
  READ(*, '(a)', IOSTAT=iostatus) f77_name
  IF (iostatus < 0) STOP               ! Halts gracefully when the names are
                                       ! read from a file and the end is reached

  IF (LEN_TRIM( f77_name ) == 0) CYCLE
  IF (INDEX(f77_name, '.') == 0) THEN
    last = LEN_TRIM(f77_name)
    f77_name(last+1:last+4) = '.for'
  END IF
  OPEN(8, file=f77_name, status='old', IOSTAT=iostatus)
  IF (iostatus /= 0) THEN
    WRITE(*, *) '** Unable to open file: ', f77_name
    CYCLE
  END IF

  pos = INDEX(f77_name, '.', BACK=.TRUE.)        ! Added BACK=.TRUE. for Unix
                                                 ! names e.g. prog.test.f
  f90_name = f77_name(1:pos) // 'f90'
  OPEN(9, file=f90_name)

!     Set up a linked list containing the lines of code

  NULLIFY(head, tail)
  ALLOCATE(head)
  tail => head
  READ(8, '(a)') head % text
  IF (head % text(1:1) == 'C' .OR. head % text(1:1) == 'c' .OR.   &
      head % text(1:1) == '*') THEN
    head % text(1:1) = '!'
  ELSE IF (head % text(1:1) == tab) THEN
    head % text = '      ' // head % text(2:)
  END IF
  head % label = ' '
  count = 1

  DO
    NULLIFY(current)
    ALLOCATE(current)
    READ(8, '(a)', IOSTAT=iostatus) current % text
    IF (iostatus /= 0) EXIT

!     Change C, c or * in column 1 to !
    IF (current % text(1:1) == 'C' .OR. current % text(1:1) == 'c' .OR.   &
        current % text(1:1) == '*'  ) THEN
      IF (LEN_TRIM(current % text) > 1) THEN
        current % text(1:1) = '!'
      ELSE
        current % text = ' '           ! Leave blank if nothing else on line
      END IF
      current % label = ' '
    ELSE
      current % label = ADJUSTL(current % text(1:5))
    END IF

    count = count + 1
    IF (current % label(1:1) == tab) THEN          ! Expand tabs
      current % label = ' '
      current % text = '      ' // current % text(2:)
    ELSE IF (current % label(1:1) == '!') THEN
      current % label = ' '
    ELSE
      current % label= ADJUSTL(current % label)
    END IF

    NULLIFY(current % next)
    tail % next => current
    tail => current
  END DO

  WRITE(*, *)'No. of lines read =', count

!---------------------------------------------------------------------------

  current => head
  NULLIFY(last_line)
  data_stmnt = .FALSE.

  DO
!     Look for blanks in columns 1-5 followed by non-blank in column 6.
!     If found, add an ampersand at the end of the previous line.

    IF (current % label == '     ' .AND. current % text(6:6) /= ' ' .AND.  &
        current % text(1:1) /= '!') THEN
      last = LEN_TRIM(last_line % text)
      last_line % text(last+3:last+3) = '&'
      current % text(6:6) = ' '
      continuation = .TRUE.
    ELSE
      data_stmnt = .FALSE.
      continuation = .FALSE.
    END IF

!     Replace tabs with single spaces
    DO
      pos = INDEX(current % text, tab)
      IF (pos == 0) EXIT
      current % text(pos:pos) = ' '
    END DO

!     Remove leading blanks
    current % text = ADJUSTL(current % text)

!     Mark regions of text which must not have their case changed.
    CALL mark_text(current % text, n_marks, pos1, pos2, continuation)

!     Convert cases of regions which are not protected.
    CALL convert_text(current % text, n_marks, pos1, pos2)

!     If line is start of a program unit, record its name
    IF (current % text(1:7) == 'PROGRAM') THEN
      prog_unit_name = current % text(1:50)
    ELSE IF (current % text(1:10) == 'SUBROUTINE') THEN
      pos = INDEX(current % text, '(') - 1
      IF (pos < 0) pos = LEN_TRIM(current % text)
      prog_unit_name = current % text(1:pos)
    ELSE IF (current % text(1:9) == 'BLOCKDATA') THEN
      prog_unit_name = current % text(1:50)
    ELSE
                               ! N.B. 'FUNCTION' could be part of a comment
      pos = INDEX(current % text, 'FUNCTION')
      IF (pos > 0 .AND. INDEX(current % text, '!') == 0 .AND.        &
                           INDEX(current % text, "'") == 0) THEN
        last = INDEX(current % text, '(') - 1
        IF (last < 0) last = LEN_TRIM(current % text)
        prog_unit_name = current % text(pos:last)
      END IF
    END IF

!     If first word is one of INTEGER, REAL, DOUBLE PRECISION, CHARACTER ,
!     LOGICAL or COMPLEX, add :: unless FUNCTION appears on the same line
!     or next non-blank character is '*' as in REAL*8.
    IF (INDEX(current % text, 'FUNCTION') == 0) THEN
      pos = 0
      IF (INDEX(current % text, 'INTEGER') == 1) THEN
        pos = 9
      ELSE IF (INDEX(current % text, 'REAL') == 1) THEN
        pos = 6
      ELSE IF (INDEX(current % text, 'DOUBLE PRECISION') == 1) THEN
        pos = 18
      ELSE IF (INDEX(current % text, 'CHARACTER') == 1) THEN
        pos = 11
      ELSE IF (INDEX(current % text, 'COMPLEX') == 1) THEN
        pos = 9
      ELSE IF (INDEX(current % text, 'LOGICAL') == 1) THEN
        pos = 9
      END IF

      IF (pos > 0) THEN
        asterisk = INDEX(current % text(pos-1:pos), '*') > 0
        IF (.NOT. asterisk) THEN
          IF (pos /= 11) THEN
            current % text = current % text(1:pos-1) // ':: ' //  &
                             ADJUSTL( current % text(pos:) )
          ELSE                         ! CHARACTER type, default length = 1
            current % text = 'CHARACTER (LEN=1) :: ' //   &
                             ADJUSTL( current % text(pos:) )
          END IF
        ELSE
          IF (pos == 11) THEN          ! CHARACTER * found
            i1 = INDEX(current % text, '*') + 1
            length = LEN_TRIM(current % text)
                                       ! Get length, could be (*)
            DO
              IF (current % text(i1:i1) /= ' ') EXIT
              IF (i1 >= length) EXIT
              i1 = i1 + 1
            END DO
            IF (current % text(i1:i1) == '(') THEN
              i1 = i1 + 1
              i2 = INDEX(current % text, ')') - 1
            ELSE
              i2 = INDEX(current % text(i1:), ' ') + i1 - 2
            END IF
            current % text = 'CHARACTER (LEN=' // current % text(i1:i2) //  &
                             ') :: ' // ADJUSTL( current % text(i2+2:) )
          END IF
        END IF
                   ! Check for 2 or more lengths in CHARACTER declaration.
                   ! e.g. CHARACTER a, b, c*10, d
                   ! Put 2nd (& later) declarations on separate lines:
                   ! CHARACTER*10 c
                   ! But check for CHARACTER*10 a(*) where last * is not a
                   ! length but a dimension
        IF (pos == 11) THEN
          pos = INDEX(current % text, '::') + 2
          DO
            i = INDEX(current % text(pos:), '*')
            IF (i == 0) EXIT
            i = i + pos - 1
            length = LEN_TRIM(current % text)
            i1 = INDEX(current % text(:i-1), ',', BACK=.TRUE.)
            i1 = MAX(pos, i1)
            i2 = INDEX(current % text(i+1:), ',')
            IF (i2 == 0) THEN
              i2 = length + 1
            ELSE
              i2 = i2 + i
            END IF
                   ! i1, i2 mark commas at beginning & end of `, name*xx,'
                   ! but we could have `name(xx, *), '
                   ! Test for * after ( , or ) before ,
            i3 = INDEX(current % text(i1:i2), '(')
            i4 = INDEX(current % text(i1:), ')')
            IF (i3 > 0)  THEN
              i4 = i4 + i1 - 1
              i2 = INDEX(current % text(i4+1:), ',')
              IF (i2 == 0) THEN
                i2 = length + 1
              ELSE
                i2 = i2 + i4
              END IF
              pos = i2 + 1
              CYCLE
            ELSE IF (i4 > 0) THEN
              i4 = i4 + i1 - 1
              i2 = INDEX(current % text(i4+1:), ',')
              IF (i2 == 0) THEN
                i2 = length + 1
              ELSE
                i2 = i2 + i4
              END IF
              pos = i2 + 1
              CYCLE
            END IF

            IF (i1 == pos .AND. i2 == length + 1) THEN
                   ! Only one declaration left on line, e.g.
                   ! CHARACTER :: name*50
              current % text = 'CHARACTER (LEN=' // current % text(i+1:length) &
                               // ') :: ' // ADJUSTL( current % text(i1:i-1) )
              EXIT
            END IF

            ALLOCATE( next_line )
            next_line % next => current % next
            current % next => next_line
            next_line % text = 'CHARACTER' // current % text(i:i2-1) //  &
                               ' ' // current % text(i1+1:i-1)
            IF (i2 < length) THEN
              current % text = current % text(:i1) // current % text(i2+1:length)
            ELSE
              current % text = current % text(:i1-1)
            END IF
          END DO
        END IF
      END IF
    END IF

!     If this is in a DATA statement, eliminate any blanks within numbers
    IF (data_stmnt .OR. current % text(1:4) == 'DATA') THEN
      CALL remove_data_blanks(current % text)
      last = LEN_TRIM(current % text)
      data_stmnt = .TRUE.
    END IF

!     If line only contains 'END', add the program unit name
    IF (LEN_TRIM(current % text) == 3 .AND. current % text(1:3) == 'END') THEN
      current % text = current % text(1:3) // ' ' // prog_unit_name
      prog_unit_name = ' '

!     Convert `enddo' to 'END DO'
    ELSE IF (current % text(1:5) == 'enddo') THEN
      current % text = 'END DO' // current % text(6:)
    END IF

    last_line => current
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

!-------------------------------------------------------------------------

!     Now convert Do-loops

  current => head
  WRITE(*, *) '      Converting DO-loops, 3-way IFs, & computed GO TOs'
  DO
    IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
      pos = INDEX(current % text, 'DO')
      IF ( pos > 0 .AND. (current % text(pos+2:pos+2) == ' ' .OR.  &
                          current % text(pos+2:pos+2) == ',' ) ) THEN
        IF ( current % text(pos+2:pos+2) == ',' )  &
                          current % text(pos+2:pos+2) = ' '
        IF (pos == 1) THEN
          OK = .TRUE.
        ELSE IF (SCAN(current % text(pos-1:pos-1), delimiters) > 0) THEN
          OK = INDEX(current % text(:pos-1), 'END ') == 0
        ELSE
          OK = .FALSE.
        END IF
        IF (OK) THEN
          text = ADJUSTL( current % text(pos+3:) )
          last = INDEX( text, ' ')
          lab = text(:last-1)
          IF (SCAN(lab(1:1), numbers) == 0) lab = ' '
          lab_length = LEN_TRIM(lab)
          IF (lab_length > 0) THEN
            pos = INDEX(lab, ',')      ! Check for a comma after label
            IF (pos > 0) THEN
               lab(pos:) = ' '
               i = INDEX(current % text, ',')
               current % text(i:i) = ' '
               lab_length = pos - 1
            END IF
            CALL do_loop_fixup(current, lab)
          END IF
        END IF

! Test for computed GO TO
      ELSE IF (INDEX(current % text, 'GO TO') > 0) THEN
        i1 = INDEX(current % text, 'GO TO')
        statement = ADJUSTL(current % text(i1+5:))
                             ! Test for a `('
        IF (statement(1:1) == '(') THEN
          OK = .TRUE.
                             ! If current line is continued, try appending
                             ! the next line
          IF (last_char(statement) == '&') THEN
            next_line => current % next
            length = LEN_TRIM(statement) + LEN_TRIM(next_line % text)
            OK = (length <= 141) .AND. (last_char(next_line % text) /= '&')
            IF (OK) THEN
              pos = LEN_TRIM(statement)
              statement = TRIM(statement(:pos-1)) // TRIM(next_line % text)
              current % next => next_line % next
              DEALLOCATE( next_line )
            END IF
          END IF

          IF (OK) THEN
                             ! Check for comma between ( and )
            pos = INDEX(statement, ')')
            IF (INDEX(statement(2:pos-1), ',') > 0) THEN
                             ! We could have something like:
                             ! IF (condition) GO TO (100, 200, 300) ivar
                             ! Before doing any more, split into 3 lines:
                             ! IF (condition) THEN
                             ! GO TO (100, 200, 300) ivar
                             ! END IF
              IF (current % text(1:2) == 'IF') THEN
                IF (current % text(3:3) == ' ' .OR.  &
                    current % text(3:3) == '(') THEN
                  current % text = current % text(:i1-1) // 'THEN'
                  i1 = 2
                  CALL insert_and_moveto_newline(current)
                  current % text = ' '
                  next_line => current
                  CALL insert_and_moveto_newline(next_line)
                  next_line % text = 'END IF'
                END IF
              END IF
                             ! Get the CASE variable or expression
              case_expr = ADJUSTL(statement(pos+1:))
              IF (case_expr(1:1) == ',') case_expr = ADJUSTL(case_expr(2:))
              current % text = current % text(:i1-1) // 'SELECT CASE ( ' // &
                               TRIM(case_expr) // ' )'
                             ! Put in pairs of lines:  CASE ( i )
                             !                         GO TO i-th label
              CALL goto_cases(statement(2:pos-1))
            END IF
          END IF
        END IF

! Look for IF, then a number as last non-blank character
      ELSE
        pos = INDEX(current % text, 'IF')
        IF (pos > 0) THEN
          last = LEN_TRIM(current % text)
          IF (SCAN(current % text(last:last), numbers) > 0) THEN
              CALL fix_3way_IF(current)
          END IF
        END IF
      END IF
    END IF

    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

!-------------------------------------------------------------------------

!     Determine INTENTs for dummy arguments

  WRITE(*, *) '      Determining INTENTs of dummy arguments'

!     Search for either FUNCTION or SUBROUTINE.
!     Extract name of program unit.

  current => head
  NULLIFY(last_line)
  outer_loop: DO
    DO
      IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
        IF (current % text(1:10) == 'SUBROUTINE') THEN
          pos = INDEX(current % text, '(') - 1
          IF (pos < 0) pos = LEN_TRIM(current % text)
          prog_unit_name = current % text(1:pos)
          EXIT
        ELSE
          pos = INDEX(current % text, 'FUNCTION')
          IF (pos > 0) THEN
            last = INDEX(current % text, '(') - 1
            IF (last < 0) last = LEN_TRIM(current % text)
            prog_unit_name = current % text(pos:last)
            EXIT
          END IF
        END IF
      END IF

      last_line => current
      current => current % next
      IF (ASSOCIATED(current, tail)) EXIT outer_loop
    END DO

!     If there is no blank line between this program unit and the previous
!     one, then insert one.

    IF (ASSOCIATED(last_line)) THEN
      IF (LEN_TRIM(last_line % text) > 0) THEN
        CALL insert_and_moveto_newline(last_line)
        last_line % text = ' '
      END IF
    END IF

    ALLOCATE( start_prog_unit )
    start_prog_unit => current

!     Find end of program unit

    DO
      current => current % next
      IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
        IF (current % text(1:3) == 'END') THEN
          IF (INDEX(current % text(5:), prog_unit_name) > 0) THEN
            ALLOCATE( end_prog_unit )
            end_prog_unit => current
            EXIT
          END IF
        END IF
      END IF
      IF (ASSOCIATED(current, tail)) EXIT outer_loop
    END DO

!     Find first & last declarations

    ALLOCATE( first_decl, last_decl )
    CALL find_declarations( start_prog_unit, end_prog_unit, first_decl, &
                            last_decl )
    IF (.NOT. ASSOCIATED(last_decl)) GO TO 100

!     Extract list of dummy arguments

    CALL get_arg_list()
    IF (numb_arg == 0) GO TO 100

!     See if the declarations contain any IMPLICIT statements

    CALL reset_defaults()
    current => first_decl
    DO
      IF( current % text(1:8) == 'IMPLICIT' ) THEN
        statement = current % text(10:)
        CALL set_implicit_types(statement)
      END IF
      IF (ASSOCIATED(current, last_decl)) EXIT
      current => current % next
    END DO

!     Search through the declarations for variable types & dimensions

    CALL get_var_types()

!     Search through rest of code to try to determine the INTENTs

    CALL get_intents()

!     Insert INTENT statements

    statement = first_decl % text
    first_decl % text = ' '
    current => first_decl
    arg => arg_start
    DO
      CALL insert_and_moveto_newline(current)
      current % text = arg % var_type
      SELECT CASE (arg % intention)
        CASE (0, 3)
          current % text = TRIM(current % text) // ', INTENT(IN OUT)'
        CASE (1)
          current % text = TRIM(current % text) // ', INTENT(IN)'
        CASE (2)
          current % text = TRIM(current % text) // ', INTENT(OUT)'
      END SELECT
      current % text = current % text(:41) // ':: ' // arg % name
      IF (arg % dim > 0)  &
          current % text = TRIM(current % text) // arg % dimensions

      IF (ASSOCIATED(arg, last_arg)) EXIT
      arg => arg % next
    END DO
    CALL insert_and_moveto_newline(current)
    current % text = statement

!     Search for, and convert, any PARAMETER statements

    current => first_decl
    DO
      IF (current % text(1:9) == 'PARAMETER') THEN
        CALL convert_parameter(current)
      END IF
      IF (ASSOCIATED(current, last_decl)) EXIT
      current => current % next
    END DO

!     Insert a blank line after the last declaration if there is not one
!     there already, or a comment.

    next_line => last_decl % next
    IF (next_line % text(1:1) /= ' ' .AND. next_line % text(1:1) /= '!') THEN
      CALL insert_and_moveto_newline(last_decl)
      last_decl % text = ' '
    END IF

!     Move onto the next SUBROUTINE or FUNCTION

    100 current => end_prog_unit
    IF (ASSOCIATED(current, tail)) EXIT
    last_line => current
    current => current % next
    IF (ASSOCIATED(current, tail)) EXIT
  END DO outer_loop

!-------------------------------------------------------------------------

!     Indenting and writing output file

!     Output header line & any continuation lines

  current => head
  continuation = .FALSE.
  DO
    IF (continuation) THEN
      WRITE(9, '(t9, a)') TRIM(current % text)
    ELSE
      WRITE(9, '(a)') TRIM(current % text)
    END IF
    ch = last_char(current % text)
    current => current % next
    IF (ch /= '&') EXIT
    continuation = .TRUE.
  END DO
!                                      Date & time stamp
  CALL DATE_AND_TIME(date, time)
  IF (ch /= ' ') WRITE(9, *)
  WRITE(9, '("! Code converted using TO_F90 by Alan Miller")')
  WRITE(9, '("! Date: ", a4, "-", a2, "-", a2, "  Time: ", a2, ":", a2,  &
        &    ":", a2)') date(1:4), date(5:6), date(7:8), time(1:2),      &
             time(3:4), time(5:6)
  IF (LEN_TRIM(current % text) > 0) WRITE(9, *)

  indent = 0
  continuation = .FALSE.
  WRITE(*, *) '      Writing file: ', f90_name

  DO
    IF (current % text(1:1) /= '!') THEN
      IF (INDEX(current % text, 'END ') > 0) THEN
        IF (INDEX(current % text, 'END SELECT') == 0) indent = MAX(indent-2, 0)
        WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
        continuation = (last_char(current % text) == '&')
      ELSE IF (INDEX(current % text, 'DO ') > 0) THEN
        WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
        continuation = (last_char(current % text) == '&')
        indent = indent + 2
                                                   ! Temporary reduction in
                                                   ! indentation for `ELSE'
      ELSE IF (INDEX(current % text, 'ELSE') > 0) THEN
        last = MAX(0, indent-2)
        WRITE(9, '(a)') blank(:last) // TRIM(current % text)
        continuation = (last_char(current % text) == '&')
                                                   ! Indent increased if `IF'
                                                   ! is followed by `THEN'
      ELSE IF (INDEX(current % text, 'IF ') > 0 .OR.           &
               INDEX(current % text, 'IF(') > 0) THEN
        current % text =  blank(:indent) // TRIM(current % text)
                                             ! If IF statement runs onto
                                             ! next line, try joining
        last = LEN_TRIM(current % text)
        IF (current % text(last:last) == '&') THEN
          next_line => current % next
          IF (last + LEN_TRIM(next_line % text) < 80) THEN
            current % text(last:last) = ' '
            current % text = TRIM(current % text) // ' ' //  &
                           TRIM(next_line % text)
            current % next => next_line % next
          END IF
        END IF

        WRITE(9, '(a)') TRIM(current % text)
        continuation = (last_char(current % text) == '&')
        next_line => current
        DO
          IF (INDEX(next_line % text, ' THEN') > 0 .OR.  &
              INDEX(next_line % text, ')THEN') > 0) THEN
            indent = indent + 2
            EXIT
          ELSE
            IF ( last_char(next_line % text) /= '&') EXIT
          END IF
          next_line => next_line % next
        END DO
      ELSE

!     If line ends with '&', attempt to join on the next line if it is short.

        last = LEN_TRIM(current % text)
        IF (last > 0) THEN
          IF (current % text(last:last) == '&') THEN
            last = LEN_TRIM(current % text(:last-1))
            next_line => current % next
            IF (last + indent + LEN_TRIM(next_line % text) < 78) THEN
              current % text = current % text(:last) // ' ' // &
                               TRIM(next_line % text)
              current % next => next_line % next
              DEALLOCATE(next_line)
            END IF
          END IF
        END IF

        IF (continuation) THEN
          WRITE(9, '(a)') blank(:indent+4) // TRIM(current % text)
        ELSE
          WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
        END IF
        continuation = (last_char(current % text) == '&')
      END IF
!     Comment line (unchanged)
    ELSE
      WRITE(9, '(a)') TRIM(current % text)
      continuation = .FALSE.
    END IF
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  CLOSE(8)
  CLOSE(9)
END DO

STOP


CONTAINS


SUBROUTINE do_loop_fixup(start, lab)

!     Convert DO-loops from:    DO xxx i=1,n    To:   DO i=1,n
!                           xxx CONTINUE              END DO

!     `start' points to the first line of the DO loop
!     `lab' is the label

TYPE (code), POINTER          :: start
CHARACTER (LEN=*), INTENT(IN) :: lab

!     Local variables

TYPE (code), POINTER :: current, end_loop
INTEGER              :: i, j, level, nmult, nl_length
LOGICAL              :: continued, jump_from_inner, referenced
CHARACTER (LEN=5)    :: label(20), next_label, text
CHARACTER (LEN=10)   :: loop_name

!-------------------------------------------------------------------
! PASS 1. Analysis
!    Find end of loop (end_loop)
!    Test for multiple loops using same label
!    Test for jumps to end of this loop from this DO loop (referenced)
!         or from inner loops (jump_from_inner)
!    Find if label is on a statement other than CONTINUE
!    Find if next executable line beyond loop is labelled (for EXIT)

current => start % next
nmult = 1
level = 0
jump_from_inner = .FALSE.
referenced = .FALSE.
DO
  IF (current % label == lab) THEN
    continued = (INDEX(current % text, 'CONTINUE') > 0)
    EXIT
  END IF

! Check for nested DO loop or multiple use of current loop

  IF (current % text(1:1) == '!' .OR. current % text(1:1) == ' ') GO TO 20
  i = INDEX(current % text, 'DO ')
  IF (i > 0 .AND. INDEX(current % text, 'END DO') == 0) THEN
    text = ADJUSTL(current % text(i+3:))
    IF (SCAN(text(1:1), numbers) > 0) THEN
      IF (text(:lab_length) == lab) THEN
        nmult = nmult + 1
      ELSE
        level = level + 1
        i = SCAN(text, ' ,')
        IF (i > 0) text = text(:i-1)
        label(level) = text
      END IF
    END IF
  END IF

! Check for end of nested loop

  IF (current % label /= '     ' .AND. level > 0) THEN
    DO
      IF (current % label == label(level)) THEN
        level = level - 1
        IF (level <= 0) EXIT
      ELSE
        EXIT
      END IF
    END DO
  END IF

! Test for GO TO current loop label

  i = INDEX(current % text, 'GO TO')
  IF (i > 0) THEN
    text = ADJUSTL(current % text(i+5:))
    IF (text(:lab_length) == lab) THEN
      IF (level > 0) THEN
        jump_from_inner = .TRUE.
      ELSE
        referenced = .TRUE.
      END IF
    END IF
  END IF

! Get next line

20 IF (.NOT. ASSOCIATED(current)) RETURN
  current => current % next
END DO

end_loop => current

! Find label of next executable line.
! First advance past any continuation lines after the end of the DO loop.

next_label = ' '
DO
  IF (last_char(current % text) /= '&') EXIT
  IF (.NOT. ASSOCIATED(current)) GO TO 10
  current => current % next
END DO

DO
  current => current % next
  IF (current % text(1:1) /= '!') EXIT
  IF (.NOT. ASSOCIATED(current)) GO TO 10
END DO
next_label = current % label
nl_length = LEN_TRIM(next_label)

!-------------------------------------------------------------------
! PASS 2. Transform beginning & end of loop

10 current => start

! Remove label from DO line
! There may be a comma after the label, if so, remove it.

i = INDEX(current % text, lab(:lab_length))
current % text = current % text(:i-1) // current % text(i+lab_length:)
length = LEN_TRIM(current % text)
DO j = i, length
  IF (current % text(j:j) == ' ') CYCLE
  IF (current % text(j:j) == ',') current % text(j:j) = ' '
  EXIT
END DO

! Jump out of inner loop detected, set up DO construct.

IF (jump_from_inner) THEN
  loop_name = 'loop' // lab
  current % text = TRIM(loop_name) // ':  ' // current % text
  current % label = ' '
END IF

! Insert `END DO' at end of loop

current => end_loop
IF (continued) THEN
  current % text = 'END DO'
  current % label = ' '
ELSE
  IF (.NOT. referenced) THEN
    current % label = ' '
    i = INDEX(current % text, lab(:lab_length))
    IF (i > 0) current % text = ADJUSTL(current % text(i+lab_length:))
  END IF
                   ! If there are continuation lines, advance to last one
  DO
    IF (last_char(current % text) == '&') THEN
      current => current % next
    ELSE
      EXIT
    END IF
  END DO
  CALL insert_and_moveto_newline(current)
  end_loop => current
  current % text = 'END DO'
END IF
IF (jump_from_inner) current % text = TRIM(current % text) // ' ' // loop_name

! Insert multiple CONTINUE's if necessary

IF (nmult > 1) THEN
  CALL insert_and_moveto_newline(current)
  end_loop => current
  current % text = lab // ' CONTINUE'
  current % label = lab
END IF

!-------------------------------------------------------------------
! PASS 3. Change GO TOs to CYCLE or EXIT where appropriate

current => start % next
IF (continued) THEN
  DO
    IF (current % text(1:1) == '!' .OR. current % text(1:1) == ' ') GO TO 30
    i = INDEX(current % text, 'GO TO')
    IF (i > 0) THEN
      text = ADJUSTL(current % text(i+5:))
      IF (text(:5) == lab) THEN
        current % text(i:) = 'CYCLE'
        IF (jump_from_inner)  &
            current % text = TRIM(current % text) // ' ' // loop_name
      ELSE IF (nl_length > 0 .AND. text(:nl_length) == next_label) THEN
        current % text(i:) = 'EXIT'
        IF (jump_from_inner)  &
            current % text = TRIM(current % text) // ' ' // loop_name
      END IF
    END IF

! Get next line

    30 current => current % next
    IF (ASSOCIATED(current, end_loop)) EXIT
    IF (.NOT.ASSOCIATED(current)) EXIT
  END DO
END IF

RETURN
END SUBROUTINE do_loop_fixup



SUBROUTINE fix_3way_IF(start)
!     Convert 3-way IFs to IF () THEN .. ELSE IF () THEN .. ELSE

TYPE (code), POINTER :: start

!     Local variables

TYPE (code), POINTER :: current
INTEGER              :: pos1, count, length, pos2, i, lab1, lab2, lab3, lenq, &
                        next_label, lenz
CHARACTER (LEN=1)    :: ch
CHARACTER (LEN=128)  :: quantity
CHARACTER (LEN=3)    :: zero_txt

current => start
length = LEN_TRIM(current % text)

!     Find closing bracket to match the opening bracket.
!     Only cases with the closing bracket on the same line are converted.

pos1 = INDEX(current % text, 'IF')

!     Check that next non-blank character after 'IF' is '('.
i = pos1 + 2
DO
  ch = current % text(i:i)
  IF (ch /= ' ') EXIT
  i = i + 1
  IF (i > length) RETURN
END DO
IF (ch /= '(') RETURN

pos1 = i
count = 1
pos2 = pos1 + 1
DO
  i = SCAN(current % text(pos2:length), '()')
  IF (i == 0) RETURN
  pos2 = i + pos2 - 1
  IF (current % text(pos2:pos2) == '(') THEN
    count = count + 1
  ELSE
    count = count - 1
  END IF
  IF (count == 0) EXIT
  pos2 = pos2 + 1
END DO

!     See if there are 3 labels after the closing bracket.

READ(current % text(pos2+1:), *, ERR=100) lab1, lab2, lab3

!     As it is probably very old code, the first alphabetic character in the
!     expression should tell us whether the quantity is REAL or INTEGER.

DO i = pos1+1, pos2-1
  ch = current % text(i:i)
  IF (ch >= 'i' .AND. ch <= 'n') THEN
    zero_txt = '0'
    lenz = 1
    EXIT
  ELSE IF (ch >= 'a' .AND. ch <= 'z') THEN
    zero_txt = '0.0'
    lenz = 3
    EXIT
  ELSE IF (i == pos2-1) THEN
    RETURN
  END IF
END DO

quantity = current % text(pos1:pos2)
lenq = LEN_TRIM(quantity)

!     Find the next executable line to see if it is labelled.
next_label = 0
DO
  IF (.NOT. ASSOCIATED(current)) EXIT
  current => current % next
  IF (current % text(1:1) == '!' .OR. LEN_TRIM(current % text) == 0) CYCLE
  IF (LEN_TRIM(current % label) > 0) READ(current % label, *) next_label
  EXIT
END DO
current => start

IF (lab1 == lab2) THEN
  current % text = current % text(:pos2-1) // ' > ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab3
  IF (lab1 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE IF (lab2 == lab3) THEN
  current % text = current % text(:pos2-1) // ' < ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  IF (lab2 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE IF (lab1 == lab3) THEN
  current % text = current % text(:pos2-1) // ' == ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  IF (lab1 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE
  current % text = current % text(:pos2-1) // ' < ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  CALL insert_and_moveto_newline(current)
  current % text = 'ELSE IF ' // quantity(1:lenq-1) // ' == ' // &
                   zero_txt(:lenz) // ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  IF (lab3 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab3
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

END IF

100 RETURN
END SUBROUTINE fix_3way_IF



SUBROUTINE insert_and_moveto_newline(current)
! Insert a new line AFTER the current line, and move `current' to point to it.

TYPE (code), POINTER :: current

!     Local variable
TYPE (code), POINTER :: new_line

ALLOCATE(new_line)
new_line % next => current % next
current % next => new_line
current => new_line

RETURN
END SUBROUTINE insert_and_moveto_newline



SUBROUTINE find_declarations( start, tail, first_decl, last_decl )
! Find the first & last declaration lines in a program unit.

TYPE (code), POINTER :: start, tail
TYPE (code), POINTER :: first_decl, last_decl

! Local variables
CHARACTER (LEN=9), PARAMETER :: declaration(13) = (/ 'IMPLICIT ', 'INTEGER  ', &
                                'REAL     ', 'DOUBLE   ', 'LOGICAL  ', &
                                'COMPLEX  ', 'DIMENSION', 'EXTERNAL ', &
                                'DATA     ', 'COMMON   ', 'PARAMETER', &
                                'SAVE     ', 'CHARACTER' /)
TYPE (code), POINTER         :: current
INTEGER                      :: pos, length, i

NULLIFY( first_decl, last_decl )

! Search for first declaration
current => start % next
search1: DO
  IF ( current % text(1:1) /= '!' .AND.  current % text(1:1) /= ' ' ) THEN
    pos = SCAN( current % text(1:13), delimiters )
    IF (pos > 0) THEN
      length = MIN(9, pos - 1)
      IF (length >= 4) THEN
        DO i = 1, 13
          IF ( current % text(:length) == declaration(i)(:length) ) THEN
            first_decl => current
            EXIT search1
          END IF
        END DO
      END IF
    END IF
  END IF

  current => current % next
  IF ( ASSOCIATED( current, tail ) ) RETURN
END DO search1

! Search for last declaration

last_decl => first_decl
DO
  IF ( current % text(1:1) /= '!' .AND.  current % text(1:1) /= ' ' ) THEN
    pos = INDEX( current % text, '=' )
    IF (pos > 0) THEN
      IF (pos < 12) RETURN
      IF (current % text(1:9) /= 'PARAMETER' .AND.  &
          current % text(1:9) /= 'CHARACTER') RETURN
    END IF

    IF ( current % text(1:4) == 'CALL' ) RETURN

    IF ( current % text(1:2) == 'IF' ) THEN
      IF ( current % text(3:3) == ' ' ) RETURN
      IF ( current % text(3:3) == '(' ) RETURN
    END IF

    IF ( current % text(1:3) == 'DO ' ) RETURN

! Skip continuation lines

    DO
      IF ( last_char(current % text) /= '&' ) EXIT
      current => current % next
    END DO

    last_decl => current
  END IF

  current => current % next
  IF ( ASSOCIATED( current, tail ) ) RETURN
END DO

RETURN
END SUBROUTINE find_declarations


SUBROUTINE get_arg_list()
! Extract list of dummy arguments

! Local variables
INTEGER :: pos, last

current => start_prog_unit
numb_arg = 0
DO                                 ! Find '(' if there are any arguments
  pos = INDEX( current % text, '(')
  IF (pos == 0) THEN
    IF ( last_char( current % text ) /= '&' ) RETURN
    current => current % next
  ELSE
    EXIT
  END IF
END DO
pos = pos + 1

NULLIFY( arg_start )
ALLOCATE( arg_start )
first_arg = .TRUE.
DO                                 ! Loop through lines of arguments
  last = SCAN(current % text(pos:), ',)')
  IF (last == 0) THEN
    IF (last_char( current % text ) /= '&' ) EXIT
    current => current % next
    pos = 1
  ELSE
    last = last + pos - 1
    NULLIFY( arg )
    ALLOCATE( arg )
    IF (first_arg) THEN
      IF (LEN_TRIM(current % text(pos:last-1)) == 0) EXIT
      arg_start => arg
      first_arg = .FALSE.
      NULLIFY( last_arg )
      ALLOCATE( last_arg )
    ELSE
      last_arg % next => arg
    END IF
    numb_arg = numb_arg + 1
    last_arg => arg

    arg % name = ADJUSTL( current % text(pos:last-1) )
    arg % intention = 0
    arg % var_type = ' '
    arg % dim = 0
    pos = last + 1
  END IF
END DO

RETURN
END SUBROUTINE get_arg_list



SUBROUTINE get_var_types()
! Search thru the declarations for the types of dummy arguments

current => first_decl
DO
  text = current % text(:30)
  IF (text(:4) == 'REAL' .OR. text(:7) == 'INTEGER' .OR.     &
      text(:6) == 'DOUBLE' .OR. text(:9) == 'CHARACTER' .OR. &
      text(:7) == 'LOGICAL' .OR. text(:7) == 'COMPLEX') THEN
                                   ! Copy the variable type to vtype
    last = INDEX(text, ' ::') - 1
    IF (last < 0) THEN
      last = INDEX(text, '*')
      IF (last == 0) THEN
        last = 24
      ELSE
        last = INDEX(text(last+2:), ' ') + last
      END IF
      i1 = last + 2
    ELSE
      i1 = last + 4
    END IF
    vtype = text(:last)
    CALL extract_declarations(i1)

  ELSE IF (text(:9) == 'DIMENSION') THEN
    i1 = 11
    vtype = ' '
    CALL extract_declarations(i1)
  END IF

  IF (ASSOCIATED(current, last_decl)) EXIT
  current => current % next
END DO

!     If there are any arguments for which the type has not been determined,
!     use the implicit types

arg => arg_start
DO
  IF (arg % var_type == ' ')   &
      arg % var_type = implicit_type(arg % name(1:1))
  IF (ASSOCIATED(arg, last_arg)) EXIT
  arg => arg % next
END DO

RETURN
END SUBROUTINE get_var_types


SUBROUTINE get_intents()
! Search thru the body of the current program unit to try to determine
! the intents of dummy arguments.

CHARACTER (LEN=80) :: last_part
INTEGER            :: j, nbrac

DO
  IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
    statement = current % text
    IF (statement(1:3) == 'IF ' .OR. statement(1:3) == 'IF(') THEN
                                       ! Split line into two parts
                                       ! IF (condition) | last_part
      i = INDEX(statement, '(')
      length = LEN_TRIM(statement)
      nbrac = 1
      DO j = i+1, length-1
        IF (statement(j:j) == ')') THEN
          nbrac = nbrac - 1
          IF (nbrac == 0) EXIT
        ELSE IF (statement(j:j) == '(') THEN
          nbrac = nbrac + 1
        END IF
      END DO
      IF (j < length) THEN
        last_part = statement(j+1:)
      ELSE
        last_part = ' '
      END IF
      statement = statement(:j)
                                       ! It is assumed that a variable inside
                                       ! an IF-expression cannot be altered
      arg => arg_start
      DO
        i = find_delimited_name(statement, arg % name)
        IF (i > 0) THEN
          IF (arg % intention == 0) arg % intention = 1
        END IF
        IF (ASSOCIATED(arg, last_arg)) EXIT
        arg => arg % next
      END DO
      statement = last_part
    END IF

    pos = INDEX(statement, '=', BACK=.TRUE.)
    IF (pos > 0) THEN
      IF (statement(pos-1:pos-1) /= '=' .AND.  &
          statement(pos-1:pos-1) /= '/' .AND.  &
          statement(pos-1:pos-1) /= '<' .AND.  &
          statement(pos-1:pos-1) /= '>') THEN

                                       ! Look for each argument name;
                                       ! is it before or after '='?
        arg => arg_start
        DO
          i = find_delimited_name(statement, arg % name)
          IF (i > 0) THEN
            IF (i < pos) THEN
              arg % intention = IOR(arg % intention, 2)
            ELSE
              IF (arg % intention == 0) arg % intention = 1
            END IF
          END IF
          IF (ASSOCIATED(arg, last_arg)) EXIT
          arg => arg % next
        END DO
      END IF
    END IF
  END IF

  IF (ASSOCIATED(current, end_prog_unit)) EXIT
  current => current % next
END DO

RETURN
END SUBROUTINE get_intents



SUBROUTINE goto_cases(text)
! Inserts pairs:
!   CASE (i)
!     GO TO i-th label
! Terminated with:
! END SELECT

CHARACTER (LEN=*), INTENT(IN OUT) :: text

INTEGER :: case_number, pos, i2

case_number = 1

DO
  pos = INDEX(text, ',')
  IF (pos > 0) THEN
    i2 = pos - 1
  ELSE
    i2 = LEN_TRIM(text)
  END IF
  CALL insert_and_moveto_newline(current)
  WRITE(current % text, '("  CASE (", i5, ")")') case_number
  CALL insert_and_moveto_newline(current)
  current % text = '    GO TO ' // TRIM(text(:i2))
  IF (pos == 0) EXIT
  text = text(pos+1:)
  case_number = case_number + 1
END DO

CALL insert_and_moveto_newline(current)
current % text = 'END SELECT'

RETURN
END SUBROUTINE goto_cases


SUBROUTINE extract_declarations(start_pos)
! Take the current line, and any continuations, look for dummy variables,
! and remove them, after storing any relevant type & dimension info.

INTEGER, INTENT(IN) :: start_pos

! Local variables

INTEGER            :: i, i1, j, ndim
CHARACTER (LEN=70) :: text

i1 = start_pos
DO
  i = SCAN(current % text(i1:), '(,')            ! Find next ( or ,
  ndim = 0
  IF (i == 0) THEN                               ! No comma or ( on this line
    IF (last_char(current % text) == '&') THEN
      current => current % next
      i1 = 1
      CYCLE
    ELSE
      text = ADJUSTL(current % text(i1:))
                                       ! Just in case there is an in-line
      pos = INDEX(text, '!')           ! comment (though illegal in F77)
      IF (pos > 0) text = text(:pos-1)

      IF (LEN_TRIM(text) == 0) RETURN
      pos = LEN_TRIM(current % text)
    END IF
  ELSE
    pos = i + i1 - 1
    IF (current % text(pos:pos) == ',') THEN     ! Comma found
      text = current % text(i1:pos-1)
    ELSE                                         ! ( found; find matching )
      count = 1
      ndim = 1
      pos = pos + 1
      DO
        j = SCAN(current % text(pos:), '(,)')
        IF (j == 0) THEN                         ! No bracket or comma
          IF (last_char(current % text) == '&') THEN
            length = LEN_TRIM(current % text)
            next_line => current % next
            current % text = TRIM(current % text(:length-1)) // ' ' //  &
                             ADJUSTL(next_line % text)
            IF (ASSOCIATED(next_line, last_decl)) last_decl => current
            current % next => next_line % next
            CYCLE
          ELSE
            RETURN
          END IF
        END IF

        pos = pos + j - 1
        SELECT CASE( current % text(pos:pos) )
          CASE ('(')
            count = count + 1
          CASE (')')
            count = count - 1
            IF (count <= 0) THEN
              text = current % text(i1:pos)
              EXIT
            END IF
          CASE (',')
            ndim = ndim + 1
        END SELECT
        pos = pos + 1
      END DO                                     ! End matching ) search
    END IF
  END IF

! Variable name isolated, with ndim dimensions
! Now see if it matches a dummy argument

  arg => arg_start
  text = ADJUSTL(text)
  IF (ndim <= 0) THEN
    length = LEN_TRIM(text)
  ELSE
    length = INDEX(text, '(') - 1
  END IF
  DO
    IF (text(:length) == arg % name) THEN        ! Argument matched
                                                 ! Insert variable type
      IF (arg % var_type == ' ') arg % var_type = vtype
      IF (ndim > arg % dim) THEN
        arg % dim = ndim
        i = INDEX(text, '(')
        arg % dimensions = text(i:)
      END IF
                                                 ! Remove variable ( & comma)
      text = ADJUSTL( current % text(pos+1:) )
      IF (LEN_TRIM(text) == 0) THEN
        IF (i1 > 1) THEN
          current % text(i1-1:) = ' '
        ELSE
          current % text = ' '
        END IF
        IF (i1 == start_pos) current % text = ' '
        RETURN
      ELSE
        IF (text(1:1) == ',') text = ADJUSTL(text(2:))
        IF (text(1:1) == '&') THEN
          next_line => current % next
          IF (i1 == start_pos) THEN
            current % text = current % text(:i1-1) // ' ' // &
                             ADJUSTL(next_line % text)
            IF (ASSOCIATED(next_line, last_decl)) last_decl => current
            current % next => next_line % next
          ELSE
            current % text = current % text(:i1-1) // '  &'
            current => next_line
            i1 = 1
          END IF
        ELSE
          current % text = current % text(:i1-1) // ' ' // text
        END IF
      END IF
      EXIT
    END IF

    IF (ASSOCIATED(arg, last_arg)) THEN
      i1 = pos + 1                               ! Skip over comma, if present
      EXIT
    END IF
    arg => arg % next
  END DO
END DO

RETURN
END SUBROUTINE extract_declarations



SUBROUTINE convert_parameter(start)

! Convert PARAMETER statements from:
! PARAMETER (name1 = value1, name2 = value2, ... )
! to:
! TYPE1, PARAMETER :: name1 = value1
! TYPE2, PARAMETER :: name2 = value2

TYPE (code), POINTER           :: start

! Local variables

TYPE (code), POINTER :: current, next_line
INTEGER              :: count, i, j, length, pos
CHARACTER (LEN=10)   :: text
CHARACTER (LEN=30)   :: vtype

current => start

! Replace opening ( with ::

i = INDEX(current % text, '(')
IF (i == 0) RETURN
current % text = TRIM(current % text(:i-1)) // ' :: ' //  &
                 ADJUSTL(current % text(i+1:))
i = INDEX(current % text, '::') + 3
DO
  j = INDEX(current % text(i:), '=')
  IF (j == 0) THEN
    IF (last_char(current % text) /= '&') RETURN
    next_line => current % next
    j = LEN_TRIM(current % text)
    current % text = TRIM(current % text(:j-1)) // next_line % text
    current % next => next_line % next
    j = INDEX(current % text(i:), '=')
    IF (j == 0) RETURN
  END IF
  j = i + j - 1
  text = ADJUSTL(current % text(i:j-1))
  CALL find_type(text, vtype, first_decl, start)

  current % text = TRIM(vtype) // ', ' // current % text
  j = j + 2 + LEN_TRIM(vtype)

! Is there another value set in this statement?
! Find end of the expression for the value, which may involve brackets
! and commas.

  10 length = LEN_TRIM(current % text)
  count = 0
  DO i = j+1,length
    SELECT CASE (current % text(i:i))
      CASE ('(')
        count = count + 1
      CASE (')')
        count = count - 1
        IF (count < 0) THEN
                             ! Remove final ) and return
          current % text = current % text(:i-1)
          RETURN
        END IF
      CASE (',')
                             ! If count = 0, there is another declaration
        IF (count == 0) THEN
                             ! Break line, check for '&' as first character
          text = ADJUSTL(current % text(i+1:))
          IF (text(1:1) == '&') THEN
            current % text = current % text(:i-1)
            current => current % next
            current % text = 'PARAMETER :: ' // ADJUSTL(current % text)
          ELSE
            ALLOCATE(next_line)
            next_line % next => current % next
            current % next => next_line
            next_line % text = 'PARAMETER :: ' // ADJUSTL(current % text(i+1:))
            current % text = current % text(:i-1)
            IF(ASSOCIATED(current, last_decl)) last_decl => next_line
            current => next_line
            start => start % next
          END IF
          EXIT
        END IF
      CASE ('&')
                             ! Expression continued on next line, merge lines
        next_line => current % next
        pos = LEN_TRIM(current % text(:i-1))
        current % text = current % text(:pos) // next_line % text
        current % next => next_line % next
        GO TO 10
    END SELECT
  END DO

  IF (i > length) EXIT
  i = 14
END DO

RETURN
END SUBROUTINE convert_parameter



SUBROUTINE find_type(vname, vtype, first_decl, last_decl)

!     Find the type of variable 'vname'

CHARACTER (LEN=*), INTENT(IN)  :: vname
CHARACTER (LEN=*), INTENT(OUT) :: vtype
TYPE (code), POINTER           :: first_decl, last_decl

! Local variables

TYPE (code), POINTER :: current
CHARACTER (LEN=30)   :: text
INTEGER              :: i1, last, length, pos

current => first_decl
length = LEN_TRIM(vname)
IF (length == 0) RETURN
DO
  text = current % text(:30)
  IF (text(:4) == 'REAL' .OR. text(:7) == 'INTEGER' .OR.     &
      text(:6) == 'DOUBLE' .OR. text(:9) == 'CHARACTER' .OR. &
      text(:7) == 'LOGICAL' .OR. text(:7) == 'COMPLEX') THEN
                                   ! Copy the variable type to vtype
    last = INDEX(text, ' ::') - 1
    IF (last < 0) THEN
      last = INDEX(text, '*')
      IF (last == 0) THEN
        last = 24
      ELSE
        last = INDEX(text(last+2:), ' ') + last
      END IF
      i1 = last + 2
    ELSE
      i1 = last + 4
    END IF
    vtype = text(:last)

! See if variable is declared on this line (& any continuation)

    DO
      pos = find_delimited_name(current % text(i1:), vname(:length))
      IF (pos == 0) THEN
        IF (last_char(current % text) == '&') THEN
          current => current % next
          i1 = 1
          CYCLE
        END IF
      END IF
      EXIT
    END DO

! Variable name found if pos > 0.

    IF (pos > 0) THEN                  ! Remove variable name
      pos = pos + i1 - 1
      current % text = current % text(:pos-1) // current % text(pos+length:)
                                       ! Delete line if only TYPE :: remains
      IF (last_char(current % text) == ':') THEN
        current % text = ' '
        RETURN
      END IF
                                       ! Remove any following comma
      i = pos
      length = LEN_TRIM(current % text)
      DO
        IF (i > length) THEN
          RETURN
        ELSE IF (current % text(i:i) == ',') THEN
          current % text = current % text(:i-1) // current % text(i+1:)
          RETURN
        ELSE IF (current % text(i:i) /= ' ') THEN
          RETURN
        END IF
        i = i + 1
      END DO
    END IF

  END IF

! If last declaration has been reached, return default type.
! Otherwise proceed to next line.

  IF (ASSOCIATED(current, last_decl)) THEN
    vtype = implicit_type(vname(1:1))
    EXIT
  ELSE
    current => current % next
  END IF
END DO

RETURN
END SUBROUTINE find_type

END PROGRAM to_f90



SUBROUTINE mark_text(text, n_marks, pos1, pos2, continuation)

!     Look for exclamation marks or quotes to find any text which must be
!     protected from case changes.
!     It is assumed that strings are NOT continued from one line to the next.
IMPLICIT NONE

CHARACTER (LEN = *), INTENT(IN)  :: text
LOGICAL, INTENT(IN)              :: continuation
INTEGER, INTENT(OUT)             :: n_marks, pos1(:), pos2(:)

!     Local variables
INTEGER                  :: mark, start, pos_exclaim, pos_sngl_quote,  &
                            pos_dbl_quote, pos, endpos
CHARACTER (LEN=1), SAVE  :: quote
LOGICAL, SAVE            :: protect = .FALSE.

mark = 1
start = 1
IF (continuation .AND. protect) THEN
  pos1(mark) = 1
  pos = 0
  GO TO 20
END IF

! Find next opening quote or exclamation mark

10 protect = .FALSE.
pos_exclaim = INDEX(text(start:80), '!')
pos_sngl_quote = INDEX(text(start:80), "'")
pos_dbl_quote = INDEX(text(start:80), '"')
IF (pos_exclaim == 0) pos_exclaim = 81
IF (pos_sngl_quote == 0) pos_sngl_quote = 81
IF (pos_dbl_quote == 0) pos_dbl_quote = 81
pos1(mark) = MIN(pos_exclaim, pos_sngl_quote, pos_dbl_quote)

IF (pos1(mark) == 81) THEN                 ! No more protected regions
  n_marks = mark - 1
  RETURN
ELSE IF (pos_exclaim == pos1(mark)) THEN   ! Rest of line is a comment
  pos1(mark) = pos1(mark) + start - 1
  pos2(mark) = 80
  n_marks = mark
  RETURN
END IF

pos = start - 1 + pos1(mark)
pos1(mark) = pos
quote = text(pos:pos)

! Search for matching quote

20 endpos = INDEX(text(pos+1:), quote)
IF (endpos > 0) THEN
  pos2(mark) = pos + endpos
  start = pos2(mark) + 1
  mark = mark + 1
  GO TO 10
END IF

! No matching end quote - it should be on the next line

pos2(mark) = 80
n_marks = mark
protect = .TRUE.

RETURN
END SUBROUTINE mark_text


SUBROUTINE convert_text(text, n_marks, pos1, pos2)

!     Convert unprotected text to upper case if it is a FORTRAN word,
!     otherwise convert to lower case.
IMPLICIT NONE

CHARACTER (LEN = *), INTENT(IN OUT) :: text
INTEGER, INTENT(IN)                 :: n_marks
INTEGER, INTENT(IN OUT)             :: pos1(:), pos2(:)

!     Local variables

INTEGER               :: length, inc = ICHAR('A') - ICHAR('a'),      &
                         pos, mark, i, i1, j, j1, j2, ptr
LOGICAL               :: matched
CHARACTER (LEN = 11)  :: fortran_word(186) = (/ "ABS        ","ACCESS     ",  &
      "ACOS       ","AIMAG      ","AINT       ","ALOG       ","ALOG10     ",  &
      "AMAX0      ","AMAX1      ","AMIN0      ","AMIN1      ","AMOD       ",  &
      "AND        ","ANINT      ","APPEND     ","ASIN       ","ASSIGN     ",  &
      "ATAN       ","ATAN2      ","BACKSPACE  ","BLANK      ","BLOCK      ",  &
      "BLOCKDATA  ","BLOCKSIZE  ","CALL       ","CCOS       ","CDABS      ",  &
      "CDCOS      ","CDEXP      ","CDLOG      ","CDSIN      ","CDSQRT     ",  &
      "CEXP       ","CHAR       ","CHARACTER  ","CLOG       ","CLOSE      ",  &
      "CMPLX      ","COMMON     ","COMPLEX    ","CONJG      ","CONTINUE   ",  &
      "COS        ","COSH       ","CSIN       ","CSQRT      ","DABS       ",  &
      "DACOS      ","DASIN      ","DATA       ","DATAN      ","DATAN2     ",  &
      "DBLE       ","DCMPLX     ","DCONJG     ","DCOS       ","DCOSH      ",  &
      "DELETE     ","DEXP       ","DIMAG      ","DINT       ","DIRECT     ",  &
      "DLOG       ","DLOG10     ","DMAX1      ","DIMENSION  ","DMIN1      ",  &
      "DMOD       ","DNINT      ","DO         ","DOUBLE     ","DSIGN      ",  &
      "DSIN       ","DSINH      ","DSQRT      ","DTAN       ","DTANH      ",  &
      "ELSE       ","ELSEIF     ","END        ","ENDFILE    ","ENDIF      ",  &
      "ENTRY      ","EQ         ","EQUIVALENCE","EQV        ","ERR        ",  &
      "EXIST      ","EXIT       ","EXP        ","EXTERNAL   ","FILE       ",  &
      "FLOAT      ","FMT        ","FORM       ","FORMAT     ","FORMATTED  ",  &
      "FUNCTION   ","GE         ","GOTO       ","GO         ","GT         ",  &
      "IABS       ","IAND       ","ICHAR      ","IDINT      ","IDNINT     ",  &
      "IEOR       ","IF         ","IFIX       ","IMPLICIT   ","INCLUDE    ",  &
      "INDEX      ","INPUT      ","INQUIRE    ","INT        ","INTEGER    ",  &
      "INTRINSIC  ","IOSTAT     ","ISIGN      ","KEEP       ","LE         ",  &
      "LEN        ","LGE        ","LGT        ","LLE        ","LLT        ",  &
      "LOG        ","LOG10      ","LOGICAL    ","LT         ","MAX        ",  &
      "MAX0       ","MAX1       ","MIN        ","MIN0       ","MIN1       ",  &
      "MOD        ","NAME       ","NAMELIST   ","NAMED      ","NE         ",  &
      "NEQV       ","NEW        ","NEXTREC    ","NONE       ","NOT        ",  &
      "NUMBER     ","OLD        ","OPEN       ","OPENED     ","OR         ",  &
      "PARAMETER  ","PAUSE      ","POSITION   ","PRECISION  ","PRINT      ",  &
      "PROGRAM    ","READ       ","REAL       ","REC        ","RECL       ",  &
      "RETURN     ","REWIND     ","SAVE       ","SCRATCH    ","SEQUENTIAL ",  &
      "SIGN       ","SIN        ","SINH       ","SNGL       ","SPACE      ",  &
      "SQRT       ","STATUS     ","STOP       ","SUBROUTINE ","TAN        ",  &
      "TANH       ","THEN       ","TO         ","TYPE       ","UNFORMATTED",  &
      "UNIT       ","UNKNOWN    ","WHILE      ","WRITE      " /)
CHARACTER (LEN = 4)   :: compare(6) = (/ ".LT.", ".LE.", ".EQ.", ".GE.",      &
                                         ".GT.", ".NE." /)
CHARACTER (LEN = 2)   :: replacement(6) = (/ "< ", "<=", "==", ">=", "> ",    &
                                             "/=" /)

!          A   B   C   D   E   F   G    H    I    J    K    L    M    N    O
!          P    Q    R    S    T    U    V    W    X    Y    Z
INTEGER, PARAMETER :: indx(27) = (/  &
           1, 20, 25, 47, 78, 92, 99, 103, 103, 121, 121, 122, 132, 139, 149, &
         153, 159, 159, 165, 177, 182, 185, 185, 187, 187, 187, 187 /)

IF (pos1(1) == 1 .AND. pos2(1) == 80) RETURN      ! Entire line protected

pos = 1
mark = 1
length = LEN_TRIM(text)
DO                                     ! Convert to upper case
  IF (n_marks >= mark .AND. pos == pos1(mark)) THEN
    pos = pos2(mark) + 1
    mark = mark + 1
    IF (pos >= length) EXIT
  END IF
  IF (text(pos:pos) >= 'a' .AND. text(pos:pos) <= 'z')           &
              text(pos:pos) = CHAR ( ICHAR(text(pos:pos)) + inc )
  pos = pos + 1
  IF (pos > length) EXIT
END DO

!     Search for `words' in text.
!     Convert to lower case if they are not FORTRAN words.
i1 = 1
pos = 1
mark = 1
DO
  IF (pos > length) EXIT
  IF (n_marks >= mark .AND. pos >= pos1(mark)) THEN
    pos = pos2(mark) + 1
    i1 = pos
    mark = mark + 1
    IF (pos >= length) EXIT
  END IF

  DO
    IF ((text(pos:pos) >= 'A' .AND. text(pos:pos) <= 'Z')        &
        .OR. (text(pos:pos) >= '0' .AND. text(pos:pos) <= '9')   &
        .OR. text(pos:pos) == '_') THEN
      pos = pos + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  pos = pos - 1
! Now i1 & pos = positions of 1st & last characters of current string

  IF (pos < i1) THEN                ! Single non-alphanumeric character
    pos = i1 + 1
    i1 = pos
    CYCLE
  END IF

  ptr = ICHAR(text(i1:i1)) - ICHAR('A') + 1
  IF (ptr < 1 .OR. ptr > 26) THEN
    pos = pos + 1
    IF (pos > length) EXIT
    i1 = pos
    CYCLE
  END IF

  matched = .FALSE.
  IF (pos > i1) THEN
    j1 = indx(ptr)
    j2 = indx(ptr+1) - 1
    DO j = j1, j2
      IF (text(i1:pos) == fortran_word(j)) THEN
        matched = .TRUE.
        EXIT
      END IF
    END DO
  END IF

! Replace .LT. with <, etc.
  IF (matched .AND. i1 > 1) THEN
    IF(text(i1-1:i1-1) == '.') THEN
      DO j = 1, 6
        IF (text(i1-1:pos+1) == compare(j)) THEN
          text(i1-1:pos+1) = ' ' // replacement(j) // ' '
          EXIT
        END IF
      END DO
      DO                                 ! Remove excess blanks
        i1 = MAX(i1, 3)
        j1 = INDEX(text(i1-2:pos+2), '  ')
        IF (j1 == 0) EXIT
        j1 = j1 + i1 - 3
        text(j1:) = text(j1+1:)
        pos2(mark) = pos2(mark) - 1      ! Adjust mark positions
        DO i = mark+1, n_marks
          pos1(i) = pos1(i) - 1
          pos2(i) = pos2(i) - 1
        END DO
        pos = pos - 1
      END DO
    END IF
  END IF

! Output line of text to screen if it contains SUBROUTINE or FUNCTION.
! Convert ENDIF to END IF, ELSEIF to ELSE IF, and GOTO to GO TO.
  IF (matched) THEN
    IF (text(i1:pos) == 'SUBROUTINE' .OR. text(i1:pos) == 'FUNCTION') THEN
      WRITE(*, '(1x, a)') text(1:length)
    ELSE IF (text(i1:pos) == 'ENDIF') THEN
      text(i1:) = 'END IF' // text(pos+1:)
      pos = pos + 1
    ELSE IF (text(i1:pos) == 'ELSEIF') THEN
      text(i1:) = 'ELSE IF' // text(pos+1:)
      pos = pos + 1
    ELSE IF (text(i1:pos) == 'GOTO') THEN
      text(i1:) = 'GO TO' // text(pos+1:)
      pos = pos + 1
    END IF
  END IF

! If text is not matched, convert to lower case, if necessary.
  IF (.NOT. matched) THEN
    DO j = i1, pos
      IF (text(j:j) >= 'A' .AND. text(j:j) <= 'Z')              &
              text(j:j) = CHAR ( ICHAR(text(j:j)) - inc )
    END DO
  END IF

  pos = pos + 1
  IF (pos > length) EXIT
  i1 = pos
END DO

RETURN
END SUBROUTINE convert_text



SUBROUTINE remove_data_blanks(text)
! Remove any blanks embedded between numerical digits in DATA statements

IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN OUT) :: text

! Local variables
INTEGER            :: length, pos, i1
CHARACTER (LEN=10) :: numbers = '1234567890'

length = LEN_TRIM(text)
i1 = 2
DO
  pos = INDEX(text(i1:length), ' ')
  IF (pos == 0) EXIT
  i1 = i1 + pos - 1
  IF (SCAN(text(i1-1:i1-1), numbers) > 0 .AND.  &
      SCAN(text(i1+1:i1+1), numbers) > 0) THEN
    text = text(:i1-1) // text(i1+1:length)
    length = length - 1
  END IF
  i1 = i1 + 2
  IF (i1 > length) EXIT
END DO

RETURN
END SUBROUTINE remove_data_blanks


FUNCTION last_char( text ) RESULT(ch)
! Return the last character on a line
IMPLICIT NONE

CHARACTER (LEN=*), INTENT(IN) :: text
CHARACTER (LEN=1)             :: ch

! Local variable
INTEGER :: last

last = LEN_TRIM( text )
IF (last == 0) THEN
  ch = ' '
ELSE
  ch = text(last:last)
END IF

RETURN
END FUNCTION last_char


FUNCTION find_delimited_name (text, name) RESULT(pos)
! Find a name in a character string with delimiters either side of it,
! or after it if it starts at position 1.
! An extended version of the intrinsic INDEX.
! pos = the position of the first character of name in text (= 0 if not found).
! N.B. When the name is short (e.g. i or n) it could occur as part of some
!      other name.

IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: text, name
INTEGER                       :: pos

! Local variables
INTEGER :: i1, ltext, lname

i1 = 1
ltext = LEN_TRIM(text)
lname = LEN_TRIM(name)
DO
  pos = INDEX(text(i1:ltext), TRIM(name))
  IF (pos == 0) RETURN
  pos = pos + i1 - 1
  IF (pos > 1) THEN
    IF ( SCAN(text(pos-1:pos-1), ' <=+-/*,') > 0 ) THEN
      IF ( SCAN(text(pos+lname:pos+lname), ' >=(+-/*,') > 0 ) RETURN
    END IF
  ELSE
    IF ( SCAN(text(pos+lname:pos+lname), ' >=(+-/*,') > 0 ) RETURN
  END IF
  i1 = pos + lname
  IF (i1 + lname > ltext) EXIT
END DO

pos = 0

RETURN
END FUNCTION find_delimited_name
