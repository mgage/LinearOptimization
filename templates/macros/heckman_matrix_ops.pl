#!/usr/bin/perl -w 

# Some WeBWorK (Perl) code for Matrix Math Objects
# by Christopher Heckman[webwork.maa.org] - Thursday, 26 December 2013, 12:05 AM
#  	I'm busy relearning how to author problems in WeBWorK (especially with the addition of MathObjects) and want to share some code that I've written. Feel free to use it; just don't claim it as your own. More will follow.
# 
# assign(M, [i, j], k) replaces M(i,j) with k.
# 
# random_vector (n, m) will generate a random vector of length n whose entries are at most m in absolute value
# 
# rand_pmDet1(n, m) will generate a random nxn matrix whose entries are bounded by a function of m (nm2, probably).
# 
# random_perm_matrixno generates a random nxn permutation matrix
# 
# zero_matrix(m,n) generates the mxn zero matrix.
# 
# rowop_multiplyrow (n, i, m)
# rowop_swaprows (n, i, j)
# rowop_addmultrow (n, i, j, m)
# 
# all generate elementary matrices (to simulate row operations); they all return nxn matrices. i and j represent relevant row numbers, and m is the multiplier.
# 
# random_rref (m, n, r, M, L)
# random_ref (m, n, r, M, L) 
# 
# generate random matrices which are in RREF or REF. The matrix returned is mxn and has rank r. The absolute values of non-pivot entries are at most M. L is a list of columns required to have pivots. Thus, 
# B = random_rref (4, 5, 3, 2, (1,3)) generates a 4x5 matrix B in RREF with rank 3, whose entries are all at most 2 in absolute value, and has pivots in columns 1 and 3. A matrix which is 4x5 and has rank 3 can be created by multiplying B on the left by an invertible 4x4 matrix (not too difficult to code; I'll probably add it later). WARNING: No common-sense checking is done, so make sure you send legitimate values (r <= min(m,n), for instance).
# 
# make_RHS (A, r, m, c?): creates a right-hand side B for a system of linear equations. A is the coefficient matrix (in REF or RREF), r is its rank, m is a bound on the size of the entries of B, and the system is made consistent iff c is nonzero. (Yes, my code for finding the rank could be used to write a version that calculates r on its own, but typically when you are writing problems, you will choose r yourself, so you may as well send along information that you already have.)
# 
# num_zeroCols (A) is the number of all-zero columns of A.
# 
# Mat2System is an item on the wish list. It takes a matrix A, a vector B, and a list of variables L, and sets up the system of linear equations nicely. Note that the order of arguments has been changed; an example of the proper way to use it is:
# 
# $s = Mat2System ($A, $b, qw(x y z));
# ...
# Consider the system of linear equations \[$s\] 
# ...
# 
# In case you don't know everything about Perl (I for one don't), the qw will take what follows and turn it into a list of strings: 
# 
# qw(x y z) ---> ("x", "y", "z");
# 
# You can even put in subscripts, as in: qw(x_1 x_2 x_3 x_4)
# 
# That's all for now.

sub assign { # assign ($L, [0,2], 13);    NOT mine
      my $self = shift; my $index = shift; my $x = shift;
      my $i = shift(@$index);
      if (scalar(@$index)) {assign($self->{data}[$i],$index,$x)}
      else {$self->{data}[$i] = Value::makeValue($x)}
    }

sub random_vector { # rand_vector (n, mult)
   my $n = shift; my $m = shift;
   my $M = Value::Matrix->I($n);
   for ($i = 0; $i < $n; $i++) { assign ($M, [$i, 0], random(-$m,$m)); }
   return Vector($M->column(1));
   }

sub rand_pmDet1 { # rand_pmDet1(n, mult)
   my $n = shift; my $m = shift;
   my $L = Value::Matrix->I($n);
   my $U = Value::Matrix->I($n);
   for ($i = 0; $i < $n; $i++) { 
      assign ($L, [$i, $i], non_zero_random(-1,1));
      for ($j = $i + 1; $j < $n; $j++) { 
         assign ($U, [$i, $j], random(-$m, $m));
         assign ($L, [$j, $i], random (-$m, $m)); 
         }
      }
   return random_perm_matrix ($n) * $L * $U;
   }

sub random_perm_matrix { # random_perm_matrixno
   my $n = shift;
   my $P = Value::Matrix->I($n);
   my $i, $j;
   my @permutation = NchooseK($n,$n);
   for ($i = 1; $i <= $n; $i++) { 
      assign($P, [$i, $i], 0); 
      assign($P, [$i, $permutation[$i]], 1);
      }
   $P;
   }

sub zero_matrix { # zero_matrix(m,n)
   my $m = shift; my $n = shift;
   my $M = Value::Matrix->I($m);
   assign ($M, [0,0], 0);
   my $A = $M -> column(1);
   for ($i = 0; $i < $m; $i++) { for ($j = 0; $j < $n; $j++) { assign ($A, [$i,$j], 0); } }
   $A
   }

 


sub rowop_multiplyrow { # rowop_multiplyrow (n, i, M);
   my $A = Value::Matrix->I(shift);
   my $i = shift;
   assign ($A, [$i - 1, $i - 1], shift);
   $A;
   }

sub rowop_swaprows { # rowop_swaprows (n, i, j);
   my $A = Value::Matrix->I(shift);
   my $i = shift; my $j = shift;
   assign ($A, [$i - 1, $i - 1], 0);
   assign ($A, [$j - 1, $j - 1], 0);
   assign ($A, [$i - 1, $j - 1], 1);
   assign ($A, [$j - 1, $i - 1], 1);
   $A;
   }

sub rowop_addmultrow { # rowop_swaprows (n, i, j, M);
   my $A = Value::Matrix->I(shift);
   my $i = shift; my $j = shift;
   my $M = shift;
   assign ($A, [$i - 1, $j - 1], $M);
   $A;
   }

sub multiple { my $m = non_zero_random(-5,5); redo if ($m == 1); $m; }


sub random_rref { # random_rref (m, n, r, M, L);
   my $m = shift; my $n = shift; my $r = shift; my $M = shift;
   my @L = @_;
   my @Lrest, $i, $indextoremove, $ro, $c;
   my $Llength = scalar @L;
   my @pivotCols;
   for ($i = 0; $i < $n; $i++) { $Lrest[$i] = $i + 1; }
   for ($i = 0; $i < $Llength; $i++) {
      $pivotCols[$i] = $L[$i];
      for ($j = 0; $j < $n - $i; $j++) 
         { if ($L[$i] == $Lrest[$j]) { $indextoremove = $j; } }
      splice (@Lrest, $indextoremove, 1);
      }
   my @restOfPivots = NchooseK ($n - $Llength, $r - $Llength);
   for ($i = 0; $i < $r - $Llength; $i++) 
      { $pivotCols[$i + $Llength] = $Lrest[$restOfPivots[$i]]; }
   @pivotCols = num_sort (@pivotCols);
   $R = Value::Matrix->I($m);
   for ($i = 0; $i < $m; $i ++) 
      { for ($j = 0; $j < $n; $j ++) { assign ($R, [$i,$j], 0); } }
   for ($i = 0; $i < $r; $i ++) { assign ($R, [$i, $pivotCols[$i] - 1], 1); }
   for ($i = 0; $i < $r - 1; $i++) {
      for ($c = $pivotCols[$i] + 1; $c < $pivotCols[$i + 1]; $c ++) { 
         for ($ro = 0; $ro <= $i; $ro++) 
            { assign ($R, [$ro, $c-1], random(-$M,$M)); } 
         }
      }
   for ($c = $pivotCols[$r - 1] + 1; $c <= $n; $c++) {
      for ($ro = 0; $ro < $r; $ro++) 
         { assign ($R, [$ro, $c - 1], random(-$M, $M)); }
      }
   $R;
   }

sub random_ref { # random_rref (m, n, r, M, L);
   my $m = shift; my $n = shift; my $r = shift; my $M = shift;
   my @L = @_;
   my @Lrest, $i, $indextoremove, $ro, $c;
   my $Llength = scalar @L;
   my @pivotCols;
   for ($i = 0; $i < $n; $i++) { $Lrest[$i] = $i + 1; }
   for ($i = 0; $i < $Llength; $i++) {
      $pivotCols[$i] = $L[$i];
      for ($j = 0; $j < $n - $i; $j++) 
         { if ($L[$i] == $Lrest[$j]) { $indextoremove = $j; } }
      splice (@Lrest, $indextoremove, 1);
      }
   my @restOfPivots = NchooseK ($n - $Llength, $r - $Llength);
   for ($i = 0; $i < $r - $Llength; $i++) 
      { $pivotCols[$i + $Llength] = $Lrest[$restOfPivots[$i]]; }
   @pivotCols = num_sort (@pivotCols);
   $R = Value::Matrix->I($m);
   for ($i = 0; $i < $m; $i ++) 
      { for ($j = 0; $j < $n; $j ++) { assign ($R, [$i,$j], 0); } }
   for ($i = 0; $i < $r; $i ++) { 
      assign ($R, [$i, $pivotCols[$i] - 1], 1); 
      for ($ro = 0; $ro < $i; $ro++) 
         { assign ($R, [$ro, $pivotCols[$i] - 1], random (-$M, $M)); }
      }
   for ($i = 0; $i < $r - 1; $i++) {
      for ($c = $pivotCols[$i] + 1; $c < $pivotCols[$i + 1]; $c ++) { 
         for ($ro = 0; $ro <= $i; $ro++) 
            { assign ($R, [$ro, $c-1], random(-$M,$M)); } 
         }
      }
   for ($c = $pivotCols[$r - 1] + 1; $c <= $n; $c++) {
      for ($ro = 0; $ro < $r; $ro++) 
         { assign ($R, [$ro, $c - 1], random(-$M, $M)); }
      }
   $R;
   }

sub make_RHS { # make_RHS (LHS, rank, multiplier, consistent?)
   my $A = shift; my $r = shift; my $m = shift; my $b = shift;
   my $mA = ($A->dimensions)[0];
   my $RHS = zero_matrix ($mA, 0);
   for ($i = 0; $i < $r; $i++) { assign ($RHS, [$i,0], random(-$m, $m)); }
   if ((! $b) && ($r < $mA)) { assign ($RHS, [$r,0], 1); }
   $RHS;
   } 

sub num_zeroCols { # num_zeroCols (A)
   my $A = shift;
   my $s, $nz, $m;
   $nz = 0;
   $m = ($A->dimensions)[0];
   for ($i = 1; $i <= ($A->dimensions)[1]; $i++) {
      $s = 0;
      for ($j = 1; $j <= $m; $j++) { $s += ($A->element($j,$i)) ** 2; }
      if ($s < 10 ** (-8)) { $nz ++; }
      } 
   $nz;
   }

sub Mat2System{ # Mat2System (A, b, qw(x y z w))
   my $coeffs = shift;
   my $vname = shift;
   my @vec = @_;
   my $srow = ($coeffs->dimensions)[0];
   my $scol = ($coeffs->dimensions)[1];
   my $vnamerow = scalar @vec;
   my $vrow = ($vname->dimensions)[0];
   die "Wrong number of rows or columns2" if  ($vrow != $srow);
   die "Wrong number of rows or columns4" if ($scol != $vnamerow);
   my $outstr="\begin{array}";
   my $s;
   $outstr .= '{r';
   for (my $j=0; $j<$scol; $j++) { $outstr .= 'rr'; }
   $outstr .= 'r}';
   for (my $j = 0; $j < $srow; $j ++) { 
      $s = 0; 
      for (my $i = 0, my $vn = 1; $i < $scol; $i++, $vn++) { 
         my $varname = $vec [$vn - 1];
         my $a=$coeffs->element($j+1,$i+1); 
         if ($a == 0) {
            if (($s > 0) || ($i < $scol - 1)) { $outstr .= '&&'; }
            else { $outstr .= '&0'; }
            }
         elsif ($a > 0) { 
            if ($a == 1) { $a = ""; } 
            if ($s==0) {$outstr .= "& $a \,$varname"; $s = 1; }
            else {$outstr .= "&+& $a \, $varname"; } 
            }
         else { 
            if ($s == 1) { 
               $a = -$a; 
               if ($a == 1) { $a = ""; } 
               $outstr .= "&- &$a \,$varname"; 
               }
            else {
               if ($a == -1) { $a = "-"; }
               $outstr = $outstr . "& $a \, $varname"; $s = 1;
               }
            } 
         } 
      $outstr .= "&=&" . $vname->element($j+1, 1). "\\"; 
      } 
   $outstr . ' \end{array}'; 
   }
1;
