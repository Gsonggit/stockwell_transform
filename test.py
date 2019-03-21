s=[]
for i in range(10000):
    tmp=input()
    if tmp == 'end':
        break
    s.append(tmp)
n=len(s[0])
m=len(s)
for i in range(m):
    for j in range(n):
        if s[i][j] == '@':
            p=i;q=j
            
def solu(s,x,y):
    if s[x][y] :
        if s[x][y] == '=':
            return solu(s,x,j+1)+solu(s,x,j-1)+solu(s,x-1,j)+solu(s,x+1,j)+1
        else :
            return 0
    else :
        return 0
print(solu(s,p,q))