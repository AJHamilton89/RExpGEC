function distance = levenshtein_distance ( s, t )

%*****************************************************************************80
%
%% LEVENSHTEIN_DISTANCE computes the Levenshtein distance between strings.
%
%  Discussion:
%
%    Let S and T be source and target strings.  Consider the task of
%    converting S to T in the minimal number of steps, involving
%    * Insertion: insert a new character;
%    * Deletion: delete a character;
%    * Substitution: swap one character for another.
%    The Levenshtein distance from S to T is the minimal number of such
%    steps required to transform S into T.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 March 2018
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    string S, T, two strings to be compared.  S is thought of as
%    the "source" string and T the "target" string.
%
%  Output:
%
%    integer DISTANCE, the Levenshtein distance between the
%    two strings.
%
  m = length ( s );
  n = length ( t );

  d = zeros ( m + 1, n + 1 );
%
%  Source prefixes can be transformed into empty string by
%  dropping all characters,
%
  for i = 1 : m
    d(i+1,1) = i;
  end
%
%  Target prefixes can be reached from empty source prefix
%  by inserting every character.
%
  for j = 1 : n
    d(1,j+1) = j;
  end

  for j = 1 : n
    for i = 1 : m
      if ( s(i) == t(j) )
        substitution_cost = 0;
      else
        substitution_cost = 1;
      end
      d(i+1,j+1) = min ( d(i,j+1) + 1, ...
                         min ( d(i+1,j) + 1, ...
                               d(i,j) + substitution_cost ) );
    end
  end
 
  distance = d(m+1,n+1);

  return
end
