from Hexago import *

#hexagon =  polygonCoordsCounterclockwise(6, (0,0), 1);
#drawBlueberry(hexagon, (0,0),16,16,4800,4800,'e',1,'H-4800-e');


p1=[(-3.0911089480061733, -20.29770503718305), (54.183730909587105, 9.846956890294585), (15.838869665994672, 12.610017391748409), (-23.00438639574441, 14.506127718598997), (13.876443158643909, 78.64804527929185), (-24.00319709782363, 14.554884066001263), (-60.88402665221193, -49.58703349469157)]
p2=[(-24.00319709782363, 14.554884066001263), (-24.003197097823634, 14.554884066001263), (-23.00438639574442, 14.506127718598998), (-23.00438639574441, 14.506127718598997)]

p3=[(0,0),(0,2)];
p4=[(-1,1),(1,1)];
print(merged(p3,p4,True,True))


'''#isSimple
p1=[(0,0),(0,1)]
p2=[(0,0),(2,0),(1,2)]
p3=[(0,0),(2,0),(2,2),(0,2)]
p4=[(0,0),(2,0),(0,2),(2,2)]
print(isSimple(p1))
print(isSimple(p2))
print(isSimple(p3))
print(isSimple(p4))'''