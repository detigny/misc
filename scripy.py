
def sinker(n):
    F = []
    for i in range(1,n+1):
        if i<10:
            f = open('/home/jetigny/JJdE/data_orange/output_0000'+str(i)+"/sink_0000"+str(i)+".csv","r+")
        else:
            f = open('/home/jetigny/JJdE/data_orange/output_000'+str(i)+"/sink_000"+str(i)+".csv","r+")
        f1 = f.readlines()
        if len(f1)==4:
            F.append([f1[2],f1[3]])
        else:
            F.append('merged')
    f.close()
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
            print(x1)
            print(x2)
            S.append(((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)**(0.5))

    return S

def correct_form(starang):
#correct for example 3E+02 to 3e2 or 1.5E-02 to 1.5e-2 without numpy
    starang = list(starang)

    k = 0
    while(True):
        if starang[k]=='E':
            break
        k = k + 1

    afterexp = starang[k:]
    afterexp[0] = 'e'
    if afterexp[1]=='+':
        afterexp.pop(1)
        if afterexp[1]=='0':
            afterexp.pop(1)
            if afterexp[1]=='0':
                afterexp = []
    else:
        if afterexp[2]=='0':
            afterexp.pop(2)
            if afterexp[2]=='0':
                afterexp = []
    return "".join(starang[:k]+afterexp)

def sink_restarter(n):

    if n<10:
        f = open('/home/jetigny/JJdE/data_orange/output_0000'+str(n)+"/sink_0000"+str(n)+".csv","r+")
    else:
        f = open('/home/jetigny/JJdE/data_orange/output_000'+str(n)+"/sink_000"+str(n)+".csv","r+")
    f1 = f.readlines()
    if len(f1)==4:
        relevant = [f1[2],f1[3]]
    else:
        relevant = 'merged'
    f.close()

    relevant = sinker(n)[n-1]
    if relevant =='merged':
        return 'The sinks are merged, do it by hand if you want to continue..'
    else:
        f= open("ic_sink_restart","w+")
        for i in range(2):
            msink = correct_form(relevant[i][16:32])
            x = str(float(correct_form(relevant[i][38:54]))-130)
            y = str(float(correct_form(relevant[i][60:76]))-130)
            z = str(float(correct_form(relevant[i][82:98]))-130)
            vx = correct_form(relevant[i][103:120])
            vy = correct_form(relevant[i][125:142])
            vz = correct_form(relevant[i][147:164])
            msmbh = correct_form(relevant[i][-28:-12])
            s = msink+','+x+','+y+','+z+','+vx+','+vy+','+vz+','+msmbh+"\n"
            f.write(s)
        f.close()
        pass
