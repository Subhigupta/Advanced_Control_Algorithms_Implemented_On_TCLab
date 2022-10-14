%========================================================================= 
% File Name     : <BPOP.m                                                  
% Usage         : BPOP_Solution = BPOP(BPOP_problem                        
% Description   : This function solves linear and quadratic mixed - integer
%                 bilevel problems of the following form                   
% min_{x1, y1}  (Q11*w1 + Q1b*w2 + c11)^T * w1 + (Q12*w2 + c12)^T * w2 +cc 
% s.t.          A11 * x1 + A12 * x2 + E11 * y1 + E12 * y2 <= b             
%               min_{x2, y2}    (Q22*w2 + Q2b*w1 + cc2)^T * w2 +(Q21*w1    
%                                                         c21)^T * w1 + cc 
%               s.t.   A21 * x1 + A22 * x2 + E21 *y1 + E22 * y2 <= b       
% w1=[x1^T y1^T]^T,   w2=[x2^T  y2^T]^                                     
% x are continuous, y are binary                                            
%------------------------------------------------------------------------- 
% Author        : Styliani Avaamidou, Efstratios N. Pistikopoulo           
% Office        : Engineering Research Building, Texas A&M University, US  
% Mail          : paroc@tamu.ed                                            
%------------------------------------------------------------------------- 
% Last Revision | Author  | Description                                     
%---------------+---------+----------------------------------------------- 
% 24-Jan-2017   | SA      | Initial Version                                 
%========================================================================= 
