
def sinker(n):
    F = []
    for i in range(1,n+1):
        if i<10:
            f = open("/home/aescala/JJdE/data/output_0000"+str(i)+"/sink_0000"+str(i)+".csv","r+")
        else:
            f = open("/home/aescala/JJdE/data/output_000"+str(i)+"/sink_000"+str(i)+".csv","r+")
        f1 = f.readlines()
        if len(f1)==4:
            F.append([f1[2],f1[3]])
        else:
            F.append('merged')
    return F

def sepper(n):
    F = sinker(n)
    S = []
    l1 = 38
    l2 = 54 
    l21 = 60
    l22 = 76
    l31 = 82
    l32 = 98

    for i in range(n):
        if F[i]=='merged':
            S.append(0)
        else:
            x1 = [float(F[i][0][l1:l2]), float(F[i][0][l21:l22]), float(F[i][0][l31:l32])]
            x2 = [float(F[i][1][l1:l2]), float(F[i][1][l21:l22]), float(F[i][1][l31:l32])]
            S.append(((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)**(0.5))

    return S
    
