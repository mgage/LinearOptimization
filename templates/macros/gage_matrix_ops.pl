#!/usr/bin/perl -w 

sub matrix_from_matrix_cols {
	my $M = shift;   # a MathObject matrix_columns
	my($n,$m) = $M->dimensions;
	my @slice = @_;
	my @columns = map {$M->column($_)->transpose->value} @slice;   
	 #create the chosen columns as rows
	 # then transform to array_refs.
	Matrix(@columns)->transpose;	#transpose and return an n by m matrix (2 dim)	
}

sub matrix_from_matrix_rows {
	my $M = shift;   # a MathObject matrix_columns
	my($n,$m) = $M->dimensions;
	my @slice = @_;
	my @rows = map {[$M->row($_)->value]} @slice;   
	 #create the chosen columns as rows
	 # then transform to array_refs.
	Matrix([@rows]); # insure that it is still an n by m matrix	(2 dim)
}

sub matrix_from_submatrix {
	my $M=shift;
	my ($row_slice,$col_slice) = @_;
	my($n,$m) = $M->dimensions;
	$row_slice=[1..$n] unless defined $row_slice;
	$col_slice=[1..$m] unless defined $col_slice;
	$M1 = matrix_from_matrix_rows($M,@$row_slice);
	return matrix_from_matrix_cols($M1, @$col_slice);
}

sub get_tableau_variable_values {
   my $mat = shift;
   my $basis =shift;
   my ($n, $m) = $mat->dimensions;
   @var = ();
   #DEBUG_MESSAGE( "start new matrix");
   foreach my $j (1..$m-2) {
      if (not $basis->contains($j)) {
            #DEBUG_MESSAGE( "j= $j not in basis");
            $var[$j-1]=0; next;
            
      } else {
            foreach my $i (1..$n-1) {
               if ( not $mat->element($i,$j) == 0 ) {
                  $var[$j-1] = ($mat->element($i,$m)/$mat->element($i,$j))->value;
                  
                  #DEBUG_MESSAGE("i=$i j=$j var = $var[$j-1] ");
                  next;
               }
             
            }
      }
  }
  push @var , ($mat->element($n,$m)/$mat->element($n,$m-1))->value;

     @var;
}

sub lp_basis_pivot {
	my ($old_tableau,$old_basis,$pivot) = @_;
	my $new_tableau= lp_clone($old_tableau);
	main::lp_pivot($new_tableau, $pivot->extract(1)-1,$pivot->extract(2)-1);	
	my $new_matrix = Matrix($new_tableau);
	my ($n,$m) = $new_matrix->dimensions;
	my $param_size = $m-$n -1;	#n=constraints+1, #m = $param_size + $constraints +2
	my $new_basis = ( $old_basis - ($pivot->extract(1)+$param_size) + ($pivot->extract(2)) )->sort;
	my @statevars = get_tableau_variable_values($new_matrix, $new_basis);
	return ( $new_tableau, $new_basis,\@statevars);
}