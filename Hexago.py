import cv2
import os
from PIL import Image, ImageDraw
import math
import copy

#variation 1 - anchor not in [coordinate, neighbor]
#variation 2 - anchor not collinear with coordinate, neighbor


def initBoard(polygon, width, height):
    im = Image.new(mode='RGBA', size=(width, height));
    #vertexCoords = polygonCoordsCounterclockwise(6,sSquare/2,sSquare/2,sHexagon); #hexagon
    drawPolygon(im, polygon, 1, 'black');
    return im;

#enclosed regions are autofilled
def makeOptimalMove(vertexCoords, variant):
    p2 = optimalMove(vertexCoords, variant);
    vertexCoords[:] = reducedApprox(merged(vertexCoords, p2, False));
    return p2;

def optimalMove(vertexCoords, variant):
    maxArea = area(vertexCoords);
    p2 = [];
    if variant=="e":
        for coordinateIndex, coordinate in enumerate(vertexCoords):
            #print(f"C: {coordinateIndex}");
            for neighborOffset in [-1, 1]:
                neighborIndex = (coordinateIndex+neighborOffset)%len(vertexCoords);
                #print(f"\tN: {neighborIndex}");
                neighbor = vertexCoords[neighborIndex];
                v = vector(neighbor, coordinate);
                mv = normalized(v);
                extendinate = summed(coordinate, mv);
                for anchorIndex, anchor in enumerate(vertexCoords):
                    if not collinearApprox(coordinate, extendinate, vertexCoords[anchorIndex]):
                        #print(f"\t\tA: {anchorIndex}");
                        debug = False;
                        if coordinateIndex==-1 and neighborIndex==0 and anchorIndex==7:
                            debug = True;
                        target = vector(anchor, summed(summed(coordinate, coordinate), mv));
                        parallelogram = [];
                        if isClockwise(coordinate, extendinate, anchor):
                            parallelogram = [anchor, extendinate, target, coordinate];
                        else:
                            parallelogram = [anchor, coordinate, target, extendinate];
                        newVertexCoords = reducedApprox(merged(vertexCoords, parallelogram, debug));
                        nArea = area(newVertexCoords);
                        if nArea > maxArea:
                            maxArea = nArea;
                            p2 = parallelogram;
    elif variant=="n":
        for coordinateIndex, coordinate in enumerate(vertexCoords):
            #print(f"C: {coordinateIndex}");
            for neighborOffset in[-1, 1]:
                neighborIndex = (coordinateIndex+neighborOffset)%len(vertexCoords);
                #print(f"\tN: {neighborIndex}");
                neighbor = vertexCoords[neighborIndex];
                for anchorIndex, anchor in enumerate(vertexCoords):
                    if not anchorIndex in [coordinateIndex, neighborIndex]:
                        #print(f"\t\tA: {anchorIndex}");
                        debug = False;
                        if coordinateIndex==-1 and neighborIndex==2 and anchorIndex==6:
                            debug = True;
                        target = vector(anchor, summed(neighbor, coordinate));
                        parallelogram = [];
                        if isClockwise(coordinate, neighbor, anchor):
                            parallelogram = [anchor, neighbor, target, coordinate];
                        else:
                            parallelogram = [anchor, coordinate, target, neighbor];
                        newVertexCoords = reducedApprox(merged(vertexCoords, parallelogram, debug));
                        nArea = area(newVertexCoords);
                        if nArea > maxArea:
                            maxArea = nArea;
                            p2 = parallelogram;
    return p2;
    
                    
def drawMove(C, E, A, vertexCoords, im):
    return;
    

def drawOptimalMoves(polygon, depth, scale, width, height, variant):
    im = initBoard(scaled(polygon, scale), width, height);
    im.show();
    for i in range(depth):
        makeOptimalMove(polygon, variant);
        drawPolygon(im, scaled(polygon, scale), 1, 'black');
        im.show();

def drawBlueberries(polygon, center, depth, scale, width, height, variant, fileName):
    centeredgon = shifted(polygon, (-center[0], -center[1]));
    centerOffset = (width/2.0, height/2.0);
    for i in range(depth):
        print(f"Drawing Blueberry {i}");
        im = Image.new(mode='RGBA', size=(width, height));
        drawPolygon(im, shifted(scaled(centergon, scale), centerOffset), 2, 'blue');
        p2 = makeOptimalMove(centeredgon, variant);
        drawPolygon(im, shifted(scaled(p2, scale), centerOffset), 2, 'green');
        absoluteDir = f"/Users/admin/Desktop/Hexago/{fileName}";
        if not os.path.exists(absoluteDir):
            os.makedirs(absoluteDir);
        im.save(f"{absoluteDir}/Blueberry {i}.png", "PNG");

def drawBlueberry(polygon, center, depth, scale, width, height, variant, fileName):
    centeredgon = shifted(polygon, (-center[0], -center[1]));
    centerOffset = (width/2.0, height/2.0);
    im = Image.new(mode='RGBA', size=(width, height));
    for i in range(depth):
        print(f"Drawing Blueberry {i}");
        drawPolygon(im, shifted(scaled(centeredgon, scale), centerOffset), 2, 'blue');
        p2 = makeOptimalMove(centeredgon, variant);
        drawPolygon(im, shifted(scaled(p2, scale), centerOffset), 2, 'green');
        absoluteDir = f"/Users/admin/Desktop/Hexago/{fileName}";
        if not os.path.exists(absoluteDir):
            os.makedirs(absoluteDir);
        im.save(f"{absoluteDir}/Blueberry Whole {i}.png", "PNG");

def centerApprox(polygon):
    if len(polygon)==0:
        return (None, None);
    xMin = polygon[0][0];
    xMax = xMin;
    yMin = polygon[0][1];
    yMax = yMin;
    for i, vertex in enumerate(polygon):
        if vertex[0]<xMin:
            xMin = vertex[0];
        elif vertex[0]>xMax:
            xMax = vertex[0];
        if vertex[1]<yMin:
            yMin = vertex[1];
        elif vertex[1]>yMax:
            yMax = vertex[1];
    return ((xMin+xMax)/2.0, (yMin+yMax)/2.0);

def betweenInclusiveSlower(A, B, C):
    if not parallel(vector(A,B), vector(B,C)):
        return False;
    if B==A or B==C or (A!=C and normalizedv(A,B)==normalizedv(A,C)):
        if magnitudev(A,B) <= magnitudev(A,C):
            return True;
    return False;

def betweenInclusiveSlow(A, B, C):
    if not collinear(A, B, C):
        return False;
    if B==A or B==C:
        return True;
    if A==C:
        return False;
    if normalizedv(A,B)==normalizedv(B,C):
        return True;
    return False;

def betweenInclusive(A, B, C):
    if not collinear(A, B, C):
        return False;
    if B[0] >= min(A[0], C[0]) and B[0] <= max(A[0], C[0]):
        if B[1] >= min(A[1], C[1]) and B[1] <= max(A[1], C[1]):
            return True;
    return False;

def betweenExclusive(A, B, C):
    if B==A or B==C:
        return False;
    return betweenInclusive(A, B, C);

def betweenInclusiveApprox(A, B, C):
    if not collinearApprox(A, B, C):
        return False;
    if areEqualApprox(B, A) or areEqualApprox(B, C):
        return True;
    if B[0] > min(A[0], C[0]) or math.isclose(B[0], min(A[0], C[0])):
        if B[0] < max(A[0], C[0]) or math.isclose(B[0], max(A[0], C[0])):
            if B[1] > min(A[1], C[1]) or math.isclose(B[1], min(A[1], C[1])):
                if B[1] < max(A[1], C[1]) or math.isclose(B[1], max(A[1], C[1])):
                    return True;
    return False;

def intersection(a1, a2, b1, b2): #if parallel, return intersection closest to a1
    if parallel(vector(a1,a2), vector(b1,b2)):
        if betweenInclusive(b1, a1, b2):
            return a1;
        if betweenInclusive(a1, b1, a2):
            if betweenInclusive(a1, b2, b1):
                return b2;
            return b1;
        if betweenInclusive(a1, b2, a2):
            if betweenInclusive(a1, b1, b2):
                return b1;
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

def intersectionApprox(a1, a2, b1, b2): #if parallel, return intersection closest to a1
    if parallelApprox(vector(a1,a2), vector(b1,b2)):
        if betweenInclusiveApprox(b1, a1, b2):
            return a1;
        if betweenInclusiveApprox(a1, b1, a2):
            if betweenInclusiveApprox(a1, b2, b1):
                return b2;
            return b1;
        if betweenInclusiveApprox(a1, b2, a2):
            if betweenInclusiveApprox(a1, b1, b2):
                return b1;
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
        if math.isclose(B, 0, abs_tol=1e-9):
            x = -C/A;
        elif math.isclose(E, 0, abs_tol=1e-9):
            x = -F/D;
        else:
            x = (F/E-C/B)/(A/B-D/E);
        if math.isclose(B, 0, abs_tol=1e-9):
            y = -D/E*x-F/E;
        else:
            y = -A/B*x-C/B;
        return (x, y);

def intersect(a1, a2, b1, b2):
    if parallel(vector(a1,a2), vector(b1,b2)):
        if betweenInclusive(a1, b1, a2):
            return True;
        if betweenInclusive(a1, b2, a2):
            return True;
        if betweenInclusive(b1, a1, b2):
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

def intersectApprox(a1, a2, b1, b2):
    if parallelApprox(vector(a1,a2), vector(b1,b2)):
        if betweenInclusiveApprox(a1, b1, a2):
            return True;
        if betweenInclusiveApprox(a1, b2, a2):
            return True;
        if betweenInclusiveApprox(b1, a1, b2):
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

def touchPolygonIndex(polygon, a1): #returns index of starting vertex of touching edge
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if betweenExclusive(vertex, a1, vnext) or vertex==a1:
            return i;
    return -1;

def touchPolygonIndexApprox(polygon, a1): #returns index of starting vertex of approximate touching edge
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if betweenInclusiveApprox(vertex, a1, vnext) and not areEqualApprox(a1, vnext):
            return i;
    return -1;
            
def intersectPolygon(polygon, a1, a2): #check if segment intersects any edge
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            return True;
    return False;

def intersectPolygonIndex(polygon, a1, a2): #returns index of starting vertex of intersecting edge
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            return i;
    return -1;

def intersectPolygonIndexNearest(polygon, a1, a2): #returns index of starting vertex of intersecting edge closest to a1
    ans = -1;
    minDistance = distance(a1, a2);
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            it = intersection(a1, a2, vertex, vnext);
            d = distance(it, a1);
            if d <= minDistance:
                minDistance = d;
                ans = i;
    return ans;

def intersectPolygonIndexNearestApprox(polygon, a1, a2): #returns index of starting vertex of approximate intersecting edge closest to a1
    ans = -1;
    minDistance = distance(a1, a2);
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersectApprox(vertex, vnext, a1, a2):
            it = intersectionApprox(a1, a2, vertex, vnext);
            d = distance(it, a1);
            if d < minDistance or math.isclose(d, minDistance, abs_tol=1e-9):
                minDistance = d;
                ans = i;
    return ans;

def intersectPolygonIndexNearestExclusive(polygon, a1, a2): #returns index of starting vertex of intersecting edge closest to but not intersecting a1
    ans = -1;
    minDistance = distance(a1, a2);
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            it = intersection(a1, a2, vertex, vnext);
            if it != a1:
                d = distance(it, a1);
                if d <= minDistance:
                    minDistance = d;
                    ans = i;
    return ans;

def intersectPolygonIndexNearestExclusiveApprox(polygon, a1, a2): #returns index of starting vertex of approximate intersecting edge closest to but not intersecting a1
    ans = -1;
    minDistance = distance(a1, a2);
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersectApprox(vertex, vnext, a1, a2):
            it = intersectionApprox(a1, a2, vertex, vnext);
            if not areEqualApprox(it, a1):
                d = distance(it, a1);
                if d < minDistance or math.isclose(d, minDistance, abs_tol=1e-9):
                    minDistance = d;
                    ans = i;
    return ans;

def intersectPolygonIndexNearestExclusives(polygon, a1, a2): #returns index of starting vertex of intersecting edge closest to a1 but not intersecting a1 or a2
    ans = -1;
    minDistance = distance(a1, a2);
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            it1 = intersection(a1, a2, vertex, vnext);
            it2 = intersection(a2, a1, vertex, vnext);
            if it1 != a1 and it2 != a2:
                d = distance(it1, a1);
                if d <= minDistance:
                    minDistance = d;
                    ans = i;
    return ans;

def intersectPolygonIndexNearestAcutest(polygon, a1, a2): #returns index of starting vertex of intersecting edge closest to and forming most acute angle with a1
    ans = -1;
    minDistance = distance(a1, a2);
    minAngle = 2*math.pi;
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        if intersect(vertex, vnext, a1, a2):
            it = intersection(vertex, vnext, a1, a2);
            d = distance(a1, it);
            if d <= minDistance:
                if d < minDistance: #closer intersection found, reset minimum angle
                    minAngle = 2*math.pi;
                ang = angle(summed(a1, vector(a2,a1)), it, vnext); #requires a0!=it, vnext!=it
                if ang <= minAngle:
                    minDistance = d;
                    minAngle = ang;
                    ans = i;
    return ans;

def angle(p1, p2, p3): #from p1 to p3, radians
    a = distance(p1, p2);
    b = distance(p3, p2);
    c = distance(p1, p3);
    #if a==0 or b==0:
    #    return -1;
    cos = (a*a+b*b-c*c)/(2*a*b);
    if math.isclose(cos, 1):
        cos = 1;
    elif math.isclose(cos, -1):
        cos = -1;
    ans = math.acos(cos);
    if isCounterclockwise(p1, p2, p3):
        return 2*math.pi-ans;
    return ans;

def longestSideLength(polygon):
    ans = 0;
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        length = distance(vertex, vnext);
        if length > ans:
            ans = length;
    return ans;

def parallel(v1, v2):
    if v1==(0,0) or v2==(0,0):
        return True;
    A = normalized(v1);
    B = normalized(v2);
    if A==B or A==(-B[0],-B[1]):
        return True;
    return False;

def parallelApprox(v1, v2):
    if areEqualApprox(v1,(0,0)) or areEqualApprox(v2,(0,0)):
        return True;
    A = normalized(v1);
    B = normalized(v2);
    if areEqualApprox(A, B) or areEqualApprox(A, (-B[0],-B[1])):
        return True;
    return False;

def collinear(A, B, C):
    v1 = vector(A, B);
    v2 = vector(B, C);
    if parallel(v1, v2):
        return True;
    return False;

def collinearApprox(A, B, C):
    if areEqualApprox(A, B) or areEqualApprox(A, C) or areEqualApprox(B, C):
        return True;
    v1 = vector(A, B);
    v2 = vector(B, C);
    if parallelApprox(v1, v2):
        return True;
    return False;
    
def polygonCoordsClockwise(n, xCenter, yCenter, radius): #fix pixel offset ========
    dtheta = 2*math.pi/n;
    theta = -math.pi/2 - dtheta/2;
    answer = [];
    for i in range(n):
        theta += dtheta;
        x = xCenter + radius*math.cos(theta);
        y = yCenter - radius*math.sin(theta);
        answer.append((x, y));
        #print(theta, answer[-1]);
    return answer;

def polygonCoordsCounterclockwise(n, xCenter, yCenter, radius):
    dtheta = 2*math.pi/n;
    theta = -math.pi/2 + dtheta*3/2;
    answer = [];
    for i in range(n):
        theta -= dtheta;
        x = xCenter + radius*math.cos(theta);
        y = yCenter - radius*math.sin(theta);
        answer.append((x, y));
    return answer;

def drawPolygon(image, vertexCoords, thickness, color):
    draw = ImageDraw.Draw(image);
    for i in range(len(vertexCoords)-1):
        draw.line(vertexCoords[i] + vertexCoords[i+1], width=thickness, fill=color);
    draw.line(vertexCoords[-1] + vertexCoords[0], width=thickness, fill=color);

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

def area(polygon): #shoelace formula
    ans = 0;
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        ans += vertex[0]*vnext[1];
        ans -= vertex[1]*vnext[0];
    ans *= 0.5;
    ans = abs(ans);
    return ans;

def merged(polygon1, polygon2, debug): #polygons valid, vertices counterclockwise
    if debug:
        print(polygon1, '\n');
        print(polygon2, '\n');
    #reduce polygons?
    #p1 = reducedApprox(polygon1);
    #p2 = reducedApprox(polygon2);
    p1 = polygon1.copy();
    p2 = polygon2.copy();
    
    #find an outer vertex as starting point
    cIndex = -1;
    nIndex = -1;
    pVertex = (0, 0);
    cVertex = (0, 0);
    nVertex = (0, 0);
    cPolygon = [];
    oPolygon = [];
    cPolygonIndex = 0;
    
    for i, vertex in enumerate(p1): #search polygon1
        if not containsInclusiveApprox(p2, vertex) or findIndexApprox(p2, vertex) != -1:
            cIndex = i;
            cVertex = vertex;
            pVertex = p1[i-1];
            cPolygon = p1;
            oPolygon = p2;
            cPolygonIndex = 1;
            break;

    if cPolygonIndex==0: #outer vertex not found in polygon1
        for i, vertex in enumerate(p2): #search polygon2
            if not containsInclusiveApprox(p1, vertex) or findIndexApprox(p1, vertex) != -1:
                cIndex = i;
                cVertex = vertex;
                pVertex = p2[i-1];
                cPolygon = p2;
                oPolygon = p1;
                cPolygonIndex = 2;
                break;
    
    #walk around current polygon, switch polygon at intersection
    ans = [cVertex];
    
    nIndex = (cIndex+1)%len(cPolygon);
    nVertex = cPolygon[nIndex];
    #print(ans[0]);

    done = False;
    while not done:
        intersectIndexExclusive = intersectPolygonIndexNearestExclusiveApprox(oPolygon, cVertex, nVertex);
        touchIndex = touchPolygonIndexApprox(oPolygon, cVertex);
        internextIndex = (intersectIndexExclusive+1)%len(oPolygon);
        touchnextIndex = (touchIndex+1)%len(oPolygon);
        #oIndex = findIndexApprox(oPolygon, cVertex);
        if debug:
            log = "";
            '''
            if oIndex != -1:
                if angle(pVertex, cVertex, nVertex) >= angle(pVertex, cVertex, oPolygon[(oIndex+1)%len(oPolygon)]):
                    log = "at vertex chan.";
                else:
                    log = "at vertex cont.";'''
            if touchIndex != -1 and not isZeroApprox(angle(pVertex, cVertex, oPolygon[touchnextIndex])) and \
            (lessExclusiveApprox(angle(pVertex, cVertex, oPolygon[touchnextIndex]), angle(pVertex, cVertex, nVertex)) \
            or (math.isclose(angle(pVertex, cVertex, oPolygon[touchnextIndex]), angle(pVertex, cVertex, nVertex), abs_tol=1e-9) and distance(cVertex, oPolygon[touchnextIndex])<distance(cVertex, nVertex))): #change polygon stay
                log = "change polygon S"
            elif intersectIndexExclusive!=-1:
                log = "change polygon M";
            else:
                log = "continue polygon";
            print(f"{f'cv: ({round(cVertex[0],2):.2f}, {round(cVertex[1],2):.2f})':<20}"\
                  f"{f'nv: ({round(nVertex[0],2):.2f}, {round(nVertex[1],2):.2f})':<20}"\
                  f"{f'pv: ({round(pVertex[0],2):.2f}, {round(pVertex[1],2):.2f})':<20}"\
                  f"{f'tnv: ({round(oPolygon[touchnextIndex][0],2):.2f}, {round(oPolygon[touchnextIndex][1],2):.2f})':<21}"\
                  f"{log:<20}cpi: {cPolygonIndex:<5}ii: {intersectIndexExclusive:<6}ci: {cIndex:<6}ti: {touchIndex:<6}");
            print("\t", cVertex);
        '''
        if oIndex != -1: #intersection at vertex
            onIndex = (oIndex+1)%len(oPolygon);
            onVertex = oPolygon[onIndex];
            if angle(pVertex, cVertex, nVertex) >= angle(pVertex, cVertex, onVertex): #change Polygon
                if areEqualApprox(onVertex, ans[0]):
                    break;
                ans.append(onVertex);
                cIndex = onIndex;
                nIndex = (cIndex+1)%len(oPolygon);
                pVertex = cVertex;
                cVertex = oPolygon[cIndex];
                nVertex = oPolygon[nIndex];
                cPolygon, oPolygon = oPolygon, cPolygon;
                cPolygonIndex = 1 if cPolygonIndex==2 else 2;
            else: #continue on cPolygon
                if areEqualApprox(nVertex, ans[0]):
                    break;
                ans.append(nVertex);
                cIndex = nIndex;
                nIndex = (cIndex+1)%len(cPolygon);
                pVertex = cVertex;
                cVertex = nVertex;
                nVertex = cPolygon[nIndex];'''
        if touchIndex != -1 and not isZeroApprox(angle(pVertex, cVertex, oPolygon[touchnextIndex])) and \
        (lessExclusiveApprox(angle(pVertex, cVertex, oPolygon[touchnextIndex]), angle(pVertex, cVertex, nVertex)) \
        or (math.isclose(angle(pVertex, cVertex, oPolygon[touchnextIndex]), angle(pVertex, cVertex, nVertex), abs_tol=1e-9) and distance(cVertex, oPolygon[touchnextIndex])<distance(cVertex, nVertex))): #change polygon stay
            tnVertex = oPolygon[touchnextIndex];
            cIndex = touchIndex;
            nIndex = touchnextIndex;
            nVertex = tnVertex;
            cPolygon, oPolygon = oPolygon, cPolygon;
            cPolygonIndex = 1 if cPolygonIndex==2 else 2;
        elif intersectIndexExclusive != -1: #change polygon move
            nVertex = intersectionApprox(cVertex, nVertex, oPolygon[intersectIndexExclusive], oPolygon[internextIndex]);
            if areEqualApprox(nVertex, ans[0]):
                break;
            ans.append(nVertex);
            if areEqualApprox(nVertex, oPolygon[internextIndex]):
                cIndex = internextIndex;
                nIndex = (cIndex+1)%len(oPolygon);
            else:
                cIndex = intersectIndexExclusive;
                nIndex = internextIndex;
            pVertex = cVertex;
            cVertex = nVertex;
            nVertex = oPolygon[nIndex];
            cPolygon, oPolygon = oPolygon, cPolygon;
            cPolygonIndex = 1 if cPolygonIndex==2 else 2;
        else: #continue on cPolygon
            if areEqualApprox(nVertex, ans[0]):
                break;
            ans.append(nVertex);
            cIndex = nIndex;
            nIndex = (cIndex+1)%len(cPolygon);
            pVertex = cVertex;
            cVertex = nVertex;
            nVertex = cPolygon[nIndex];
    return ans;

def reduced(polygon): #polygon valid (no self-intersections); combines collinear vertices, removes duplicate vertices
    ans = [];
    for i, cv in enumerate(polygon):
        pv = ans[-1] if len(ans) else polygon[i-1];
        nv = polygon[(i+1)%len(polygon)];
        if betweenExclusive(pv, cv, nv) or cv==pv:
            continue;
        ans.append(cv);
    return ans;

def reducedApprox(polygon): #polygon valid (no self-intersections); combines collinearApprox vertices, removes areEqualApprox vertices
    ans = [];
    for i, cv in enumerate(polygon):
        pv = ans[-1] if len(ans) else polygon[i-1];
        nv = polygon[(i+1)%len(polygon)];
        if areEqualApprox(cv, pv): #don't include cv
            continue;
        if betweenInclusiveApprox(pv, cv, nv) and not areEqualApprox(cv, nv):
            continue;
        ans.append(cv);
    return ans;
        
def findIndex(array, value):
    for i, v in enumerate(array):
        if v==value:
            return i;
    return -1;

def findIndexApprox(polygon, vertex):
    for i, v in enumerate(polygon):
        if areEqualApprox(v, vertex):
            return i;
    return -1;

def areEqualApprox(v1, v2):
    if math.isclose(v1[0], v2[0], abs_tol=1e-9) and math.isclose(v1[1], v2[1], abs_tol=1e-9):
        return True;
    return False;

def isZeroApprox(n):
    if math.isclose(n, 0, abs_tol=1e-9):
        return True;
    return False;

def swap(polygon1, polygon2): #doesn't work
    temp = polygon1.copy();
    polygon1 = polygon2.copy();
    polygon2 = temp;

def lessExclusiveApprox(a, b):
    if a<b and not math.isclose(a, b, abs_tol=1e-9):
        return True;
    return False;

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

def containsExclusiveApprox(polygon, point): #polygon vertices counterclockwise
    if not containsConvexExclusiveApprox(polygon, point):
        return False;

    #cut all concave vertices, check if each triangle contains point
    vcoords = polygon.copy();
    done = False;
    while not done: #vcoords not convex
        done = True;
        for i, vertex in enumerate(vcoords):
            vprev = vcoords[i-1];
            vnext = vcoords[(i+1)%len(polygon)];
            if isClockwise(vprev, vertex, vnext): #concavity at vertex
                triangle = [vprev, vnext, vertex];
                if containsConvexInclusiveApprox(triangle, point):
                    return False;
                #remove concavity
                vcoords.pop(i);
                done = False;
                break;
    #vcoords now convex
    if containsConvexExclusiveApprox(vcoords, point):
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

def containsInclusiveApprox(polygon, point): #polygon vertices counterclockwise
    if not containsConvexInclusiveApprox(convexed(polygon), point):
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
                if containsConvexInclusiveApprox(triangle, point) and not betweenInclusiveApprox(vprev, point, vertex) and not betweenInclusiveApprox(vertex, point, vnext):
                    return False;
                #point not in concavity, cut concavity away
                vcoords.pop(i);
                done = False;
                break;
    #vcoords now convex
    if containsConvexInclusiveApprox(vcoords, point):
        return True;
    return False;

def containsConvexExclusive(polygon, point): #polygon convex, vertices counterclockwise
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        if not isCounterclockwise(vertex, vnext, point):
            return False;
    return True;

def containsConvexExclusiveApprox(polygon, point): #all points in polygon not on or close to boundary
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        if isClockwiseApprox(vertex, vnext, point):
            return False;
    return True;

def containsConvexInclusive(polygon, point): #polygon convex, vertices counterclockwise
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext = polygon[(i+1)%len(polygon)];
        if isClockwise(vertex, vnext, point):
            return False;
    return True;

def containsConvexInclusiveApprox(polygon, point): #polygon convex, vertices counterclockwise
    for i in range(len(polygon)):
        vertex = polygon[i];
        vnext =  polygon[(i+1)%len(polygon)];
        if not isCounterclockwiseApprox(vertex, vnext, point):
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

def isCounterclockwiseApprox(A, B, C):
    BA = vector(B, A);
    BC = vector(B, C);
    if parallelApprox(BA, BC):
        return True;
    argBA = math.atan2(BA[1], BA[0]);
    argBC = math.atan2(BC[1], BC[0]);
    dth = (argBC-argBA+2*math.pi)%(2*math.pi);
    if dth > math.pi or math.isclose(dth, math.pi):
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

def isClockwiseApprox(A, B, C):
    BA = vector(B, A);
    BC = vector(B, C);
    if parallelApprox(BA, BC):
        return True;
    argBA = math.atan2(BA[1], BA[0]);
    argBC = math.atan2(BC[1], BC[0]);
    dth = (argBC-argBA+2*math.pi)%(2*math.pi);
    if dth < math.pi or math.isclose(dth, math.pi):
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

def distance(A, B):
    return magnitude(vector(A, B));

def summed(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1]);

def perimeter(polygon):
    ans = 0;
    for i, vertex in enumerate(polygon):
        vnext = polygon[(i+1)%len(polygon)];
        ans += distance(vertex, vnext);
    return ans;

def rounded(tup, digits):
    ans = tuple([round(i, digits) for i in tup]);
    return ans;
