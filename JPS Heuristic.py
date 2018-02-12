from random import *
from operator import *
from math import *
seed(a=10000)
def individual_x1(mi,ma):
    x1=uniform(mi,ma)
    return x1

def individual_x2(mi,ma):
    x2=uniform(mi,ma)
    return x2

def individual_x3(mi,ma):
    x3=uniform(mi,ma)
    return x3

def individual_x4(mi,ma):
    x4=uniform(mi,ma)
    return x4

def individual_lambda(mi,ma):
    lambd=uniform(mi,ma)
    return lambd

def individual_rc(mi,ma):
    rc=uniform(mi,ma)
    return rc

def individual_b(mi,ma):
    b=uniform(mi,ma)
    return b

def fitness(x1,x2,x3,x4,lambd,rc,b):
    t,rha,rhc,pa,pc,a,ilimitden,l,ns=353.15,1,1,3,5,27,0.860,0.127,24
    sum=0
    for x in range(0,15):
        i=1.1+x
        iden=i/a
        psath2o=10**((0.0295*(t-273.15))-(0.0000919*((t-273.15)**2))+(0.000000144*((t-273.15)**3))-2.18)
        ph2=0.5*rha*psath2o*((1/((rha*psath2o/pa)*exp((1.635*iden)/t**1.334)))-1)
        #print(ph2)
        po2=rhc*psath2o*((1/((rhc*psath2o/pc)*exp((4.192*iden)/t**1.334)))-1)
        #print(po2)
        enernst=1.229-(0.85*0.001*(t-298.15))+(0.000043085*t*log1p(ph2*sqrt(po2)))
        #print(enernst)
        co2=po2/(5080000*exp(-498/t))
        nact_model=-(x1+(0.003*t)+(x3*t*log1p(co2))+(x4*t*log1p(i)))
        #print(x1,(0.003*t),x3*t*log1p(co2),x4*t*log1p(i))
        nconc_model=-0.035*log1p(1-(iden/ilimitden))
        #print(b)
        pm_model=(181.6*(1+(0.03*iden)+0.062*(t/303)*(iden**2.5)))/((lambd-0.634-3*iden)*exp(4.18*((t-303)/t)))
        rm_model=(pm_model*l)/a
        nohm_model=i*(rm_model+rc)
        vcell_model=enernst-nact_model-nohm_model-nconc_model
        #print(enernst)
        #print(nact_model)
        #print(nohm_model)
        #print(nconc_model)
        vstack_model=ns*vcell_model
        #print(vstack_model)             
        nact_actual=-(-0.944957+(0.00301801*t)+(0.00007401*t*log1p(co2))+(-0.000188*t*log1p(i)))
        #print(-0.944957,0.00301801*t,0.00007401*t*log1p(co2),-0.000188*t*log1p(i))
        nconc_actual=-0.02914489*log1p(1-(iden/ilimitden))
        pm_actual=(181.6*(1+(0.03*iden)+0.062*(t/303)*(iden**2.5)))/((23-0.634-3*iden)*exp(4.18*((t-303)/t)))
        rm_actual=(pm_model*l)/a
        nohm_actual=i*(rm_model+0.0001)
        vcell_actual=enernst-nact_actual-nohm_actual-nconc_actual
        vstack_actual=ns*vcell_actual
        #print(enernst,nact_actual,nohm_actual,nconc_actual)
        #print(enernst,nact_model,nohm_model,nconc_model)
        #print(vstack_actual)
        sum+=(vstack_actual-vstack_model)**2
    #print(sum)
    return sum


def pop(itr=100):
    lst=[]
    for x in range(itr):
        innerlst=[]
        x1=individual_x1(-0.952,-0.944)
        x2=individual_x2(0.001,0.005)
        x3=individual_x3(0.000074,0.000078)
        x4=individual_x4(-0.000198,-0.0000188)
        lambd=individual_lambda(14,23)
        rc=individual_rc(0.0001,0.0008)
        b=individual_b(0.016,0.5)
        innerlst.extend([x1,x2,x3,x4,lambd,rc,b])
        lst.append(innerlst)
    return lst


def evals(e):
    return e-1

def avg_fitness():
    lst=pop()
    #print(parents)
    #print(len(lst))
    sum_fitness=0
    for x in range(len(lst)):
        x1=lst[x][0]
        x2=lst[x][1]
        x3=lst[x][2]
        x4=lst[x][3]
        lambd=lst[x][4]
        rc=lst[x][5]
        b=lst[x][6]
        sum_fitness+=fitness(x1,x2,x3,x4,lambd,rc,b)
    return (sum_fitness/(len(lst)))

def mutation(a,mutation=0.7):
    x=0.3
    old_val=0
    new_val=0
    hi=0
    lo=0
    v=0
    if mutation>random():
        pos_to_mutate=randint(0,len(a)-1)
        if pos_to_mutate==0:
            old_val,hi,lo=individual_x1(-0.952,-0.944),-0.944,-0.952
        if pos_to_mutate==1:
            old_val,hi,lo=individual_x2(0.001,0.005),0.005,0.001
        if pos_to_mutate==2:
            old_val,hi,lo=individual_x3(0.000074,0.000078),0.000078,0.000074
        if pos_to_mutate==3:
            old_val,hi,lo=individual_x4(-0.000198,-0.0000188),-0.0000188,-0.000198
        if pos_to_mutate==4:
            old_val,hi,lo=individual_lambda(14,23),23,14
        if pos_to_mutate==5:
            old_val,hi,lo=individual_rc(0.0001,0.0008),0.0008,0.0001
        if pos_to_mutate==6:
            old_val,hi,lo=individual_b(0.016,0.5),0.5,0.016

        v=x*(hi-lo)*random()
        if random()>0.5:
            v=-v
        else:
            v=-v
        new_val=old_val+v
        if new_val<lo:
            if old_val==lo:
                new_val=lo
            else:
                new_val=lo+uniform(0,1)*(old_val-lo)
        elif new_val>hi:
            if old_val==hi:
                new_val=hi
            else:
                new_val=hi-uniform(0,1)*(hi-old_val)
        a[pos_to_mutate]=new_val
    return a
        
            
            

def jps():
    lst=pop()
    davg=avg_fitness()
    a1=choice(lst)
    e=4000
    k=1

    a=a1
    a_best=a1

    while e>0:
        a=move(a,lst,e,davg)
        if individual_fitness(a)<individual_fitness(a1):
            a_best=a
        e=evals(e)
        k+=1
    return a_best

def move(a,lst,e,davg):
    b=mutation(a)
    d=individual_fitness(b)-individual_fitness(a)
    if d<0:
        return b
    else:
        p=((0.99998)**e)/(1+(d/davg)**2)
        return b if random()<p else a


def individual_fitness(a):
    x1=a[0]
    x2=a[1]
    x3=a[2]
    x4=a[3]
    lambd=a[4]
    rc=a[5]
    b=a[6]
    return(fitness(x1,x2,x3,x4,lambd,rc,b))


def main():
    #seed(a=2)
    a_best=jps()
    fitness_value=individual_fitness(a_best)
    print(fitness_value)

main()
