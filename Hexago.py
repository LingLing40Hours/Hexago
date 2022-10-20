import cv2
import os
from PIL import Image, ImageDraw
import math
import copy


def initBoard(sHexagon, sSquare):
    im = Image.new(mode='RGBA', size=(sSquare,sSquare));
    vertexCoords = polygonCoords(6,sSquare/2,sSquare/2,sHexagon);
    drawPolygon(im, vertexCoords, 1);
    return im;

#enclosed regions are autofilled
def optimalMove(vertexCoords):
    for coordinateIndex in range(len(vertexCoords)):
        coordinate = vertexCoords[coordinateIndex];
        for neighborOffset in range(-1, 2, 2):
            neighborIndex = (coordinateIndex+neighborOffset)%len(vertexCoords);
            neighbor = vertexCoords[neighborIndex];
            v = coordinate-neighbor;
            mv = normalized(v);
            extendinate = coordinate+v;
            for anchorIndex in range(len(vertexCoords)):
                if anchorIndex!=coordinateIndex and anchorIndex!=neighborIndex:
                    target = 2*coordinate+mv-anchor;
                    #target inside polygon
                    #target outside polygon
                    
def drawMove(C, E, A, vertexCoords, im):
    return;

def intersection(a1, a2, b1, b2):
    if parallel(vector(a1,a2), vector(b1,b2)):
        if normalizedv(a1,b1)==normalizedv(a1,a2) and magnitudev(a1,b1)<=magnitudev(a1,a2): #b1 between a1, a2 inclusive
            return b1;
        if normalizedv(a1,b2)==normalizedv(a1,a2) and magnitudev(a1,b2)<=magnitudev(a1,a2): #b2 between a1, a2 inclusive
            return b2;
        return (None, None);
    else: #Ax+By+C=0, Dx+Ey+F=0
        A = a1[1]-a2[1];
        B = a2[0]-a1[0];
        C = -B*a1[1]-A*a1[0];
        D = b1[1]-b2[1];
        E = b2[0]-b1[0];
        F = -E*b1[1]-D*b1[0];
        x=0;
        y=0;
        if B==0:
            x = -C/A;
        elif E==0:
            x = -F/D;
        else:
            x = (F/E-C/B)/(A/B-D/E);
        if B==0:
            y = -D/E*x-F/E;
        else:
            y = -A/B*x-C/B;
        return (x, y);
        

def intersect(a1, a2, b1, b2):
    if parallel(vector(a1,a2), vector(b1,b2)):
        if normalizedv(a1,b1)==normalizedv(a1,a2) and magnitudev(a1,b1)<=magnitudev(a1,a2): #b1 between a1, a2 inclusive
            return True;
        if normalizedv(a1,b2)==normalizedv(a1,a2) and magnitudev(a1,b2)<=magnitudev(a1,a2): #b2 between a1, a2 inclusive
            return True;
        return False;
    cha1a2b1 = chirality(a1,a2,b1);
    cha1a2b2 = chirality(a1,a2,b2);
    cha2a1b1 = chirality(a2,a1,b1);
    cha2a1b2 = chirality(a2,a1,b2);
    chb1b2a1 = chirality(b1,b2,a1);
    chb1b2a2 = chirality(b1,b2,a2);
    chb2b1a1 = chirality(b2,b1,a1);
    chb2b1a2 = chirality(b2,b1,a2);
    if cha1a2b1!=-1 and cha1a2b2!=1 and chb1b2a2!=-1 and chb1b2a1!=1:
        return True;
    if cha1a2b2!=-1 and cha1a2b1!=1 and chb1b2a1!=-1 and chb1b2a2!=1:
        return True;
    return False;
            
        

def parallel(v1, v2):
    if v1==(0,0) or v2==(0,0):
        return True;
    A = normalized(v1);
    B = normalized(v2);
    if A==B or A==(-B[0],-B[1]):
        return True;
    return False;

def collinear(A, B, C):
    v1 = vector(A, B);
    v2 = vector(B, C);
    if parallel(v1, v2):
        return True;
    return False;
            
    
def polygonCoords(n, xCenter, yCenter, radius): #fix pixel offset ========
    dtheta = 2*math.pi/n;
    theta = -math.pi/2 - dtheta/2;
    answer = [];
    for i in range(n):
        theta += dtheta;
        x = xCenter + radius*math.cos(theta);
        y = yCenter - radius*math.sin(theta);
        answer.append((x, y));
        print(theta, answer[-1]);
    return answer;

def drawPolygon(image, vertexCoords, thickness):
    draw = ImageDraw.Draw(image);
    for i in range(len(vertexCoords)-1):
        draw.line(vertexCoords[i] + vertexCoords[i+1], width=thickness, fill='black');
    draw.line(vertexCoords[-1] + vertexCoords[0], width=thickness, fill='black');

def drawCircle(image, point, r, borderColor, borderWidth):
    draw = ImageDraw.Draw(image);
    box = (point[0]-r, point[1]-r, point[0]+r, point[1]+r);
    draw.ellipse(box, outline=borderColor, width=borderWidth, fill=(0,0,0,0));

def drawDot(image, point, r, RGBA):
    draw = ImageDraw.Draw(image);
    box = (point[0]-r, point[1]-r, point[0]+r, point[1]+r);
    draw.ellipse(box, outline=RGBA, width=1, fill=RGBA);

def magnitude(v):
    a, b = v;
    return math.sqrt(a*a+b*b);

def normalized(v):
    m = magnitude(v);
    a, b = v;
    #if m==0:
    #    return (0, 0);
    return (a/m, b/m);

def normalizedv(a, b):
    return normalized(vector(a, b));

def magnitudev(a, b):
    return magnitude(vector(a, b));

def area(polygon):
    ans = 0;
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        ans += vertex[0]*vnext[1];
        ans -= vertex[1]*vnext[0];
    ans *= 0.5;
    ans = abs(ans);
    return ans;

def merged(polygon1, polygon2): #polygons valid and counterclockwise
    #find an outer vertex as starting point
    start = (0, 0);
    startIndex = 0;
    startPolygon = 0;
    for i, vertex in enumerate(polygon1):
        if not contains(polygon2, vertex):
            start = vertex;
            startIndex = i;
            startPolygon = 1;
            break;
    if startPolygon==0: #polygon2 contains polygon1
        for i, vertex in enumerate(polygon2):
            if not contains(polygon1, vertex):
                start = vertex;
                startIndex = i;
                startPolygon = 2;
                break;

    #walk around polygon while checking for intersection
    

def convexed(polygon): #vertices counterclockwise
    vcoords = polygon.copy();
    done = False;
    while not done:
        done = True;
        for i, vertex in enumerate(vcoords):
            vprev = vcoords[i-1];
            vnext = vcoords[(i+1)%len(vcoords)];
            if isClockwise(vprev, vertex, vnext):
                vcoords.pop(i);
                done = False;
                break;
    return vcoords;

def convexedStrict(polygon): #vertices counterclockwise; result interior angles less than 180
    vcoords = polygon.copy();
    done = False;
    while not done:
        done = True;
        for i, vertex in enumerate(vcoords):
            vprev = vcoords[i-1];
            vnext = vcoords[(i+1)%len(vcoords)];
            if not isCounterclockwise(vprev, vertex, vnext):
                vcoords.pop(i);
                done = False;
                break;
    return vcoords;

def containsExclusive(polygon, point): #polygon vertices counterclockwise
    if not containsConvexExclusive(convexed(polygon), point):
        return False;
    
    #cut all concave vertices, check if each triangle contains point
    vcoords = polygon.copy();
    done = False;
    while not done:
        done = True;
        for i, vertex in enumerate(vcoords):
            vprev = vcoords[i-1];
            vnext = vcoords[(i+1)%len(vcoords)];
            if isClockwise(vprev, vertex, vnext): #concavity at vertex
                triangle = [vprev, vnext, vertex];
                if containsConvexInclusive(triangle, point):
                    return False;
                #cut vertex away to extend vcoords
                vcoords.pop(i);
                done = False;
                break;
    #vcoords now convex
    if containsConvexExclusive(vcoords, point):
        return True;
    return False;

def containsInclusive(polygon, point): #polygon vertices counterclockwise
    if not containsConvexInclusive(convexed(polygon), point):
        return False;

    #cut all concave vertices, check if each triangle contains point
    vcoords = polygon.copy();
    done = False;
    while not done:
        done = True;
        for i, vertex in enumerate(vcoords):
            vprev = vcoords[i-1];
            vnext = vcoords[(i+1)%len(vcoords)];
            if isClockwise(vprev, vertex, vnext): #concavity at vertex
                triangle = [vprev, vnext, vertex];
                if containsConvexExclusive(triangle, point):
                    return False;
                if containsConvexInclusive(triangle, point) and collinear(vprev, point, vnext) and point!=vnext and point!=vprev:
                    return False;
                #cut vertex away to extend vcoords
                vcoords.pop(i);
                done = False;
                break;
    #vcoords now convex
    if containsConvexInclusive(vcoords, point):
        return True;
    return False;

def containsConvexExclusive(polygon, point): #polygon convex, vertices counterclockwise
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        if not isCounterclockwise(vertex, vnext, point):
            return False;
    return True;

def containsConvexInclusive(polygon, point): #polygon convex, vertices counterclockwise
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        if isClockwise(vertex, vnext, point):
            return False;
    return True;

def dot(v1, v2):
    return v1[0]*v2[0]+v1[1]*v2[1];

def isCounterclockwise(A, B, C):
    BA = vector(B, A);
    BC = vector(B, C);
    if parallel(BA, BC):
        return False;
    argBA = math.atan2(BA[1], BA[0]);
    argBC = math.atan2(BC[1], BC[0]);
    dth = (argBC-argBA+2*math.pi)%(2*math.pi);
    if dth > math.pi:
        return True;
    return False;
    
def isClockwise(A, B, C):
    BA = vector(B, A);
    BC = vector(B, C);
    if parallel(BA, BC):
        return False;
    argBA = math.atan2(BA[1], BA[0]);
    argBC = math.atan2(BC[1], BC[0]);
    dth = (argBC-argBA+2*math.pi)%(2*math.pi);
    if dth < math.pi:
        return True;
    return False;

def chirality(A, B, C): #clockwise 1, counterclockwise -1, collinear 0
    if isClockwise(A, B, C):
        return 1;
    if isCounterclockwise(A, B, C):
        return -1;
    return 0;
    

def vector(A, B):
    return (B[0]-A[0], B[1]-A[1]);

def shifted(polygon, v):
    ans = polygon.copy();
    for i in range(len(ans)):
        x, y = ans[i];
        ans[i] = (x+v[0], y+v[1]);
    return ans;

def scaled(polygon, k):
    ans = polygon.copy();
    for i in range(len(ans)):
        x, y = ans[i];
        ans[i] = (k*x, k*y);
    return ans;
