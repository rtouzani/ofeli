<?xml version="1.0" encoding="ISO-8859-1" ?>
<OFELI_File>
<info>
   <title>Finite Element Mesh of a rectangle</title>
   <date>January 1, 2008</date>
   <author>R. Touzani</author>
</info>
<Mesh dim="2">
   <Nodes>
     0.00   1.00    1        0.00   0.75    0        0.00   0.50    0        0.00   0.25    0
     0.00   0.00    1        1.00   1.00    1        1.00   0.75    0        1.00   0.50    0
     1.00   0.25    0        1.00   0.00    1        2.00   1.00    1        2.00   0.75    0
     2.00   0.50    0        2.00   0.25    0        2.00   0.00    1        3.00   1.00    1
     3.00   0.75    0        3.00   0.50    0        3.00   0.25    0        3.00   0.00    1
   </Nodes>
   <Elements shape="triangle" nodes="3">
        1    2    6    1        2    7    6    1        2    3    7    1        3    8    7    1
        3    4    8    1        4    9    8    1        4    5    9    1        5   10    9    1
        6    7   11    1        7   12   11    1        7    8   12    1        8   13   12    1
        8    9   13    1        9   14   13    1        9   10   14    1       10   15   14    1
       11   12   16    1       12   17   16    1       12   13   17    1       13   18   17    1
       13   14   18    1       14   19   18    1       14   15   19    1       15   20   19    1
    </Elements>
</Mesh>
</OFELI_File>

