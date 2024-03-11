import random as r
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe

etp1 =  [[('red',[24000,24000],[57000,67000])],[('red',[27000,27000],[57000,67000])]]
etp2 = [[('blue',[33000,33000],[57000,67000])],[('blue',[36000,36000],[57000,67000])]]
stp1,stp2=etp1,etp2

def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    plt.grid()
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(tp1)
    draw_pair(tp2)
    plt.text(30000-1000,61500, 'X', size=20)
    return



#RILs and NILs

etp1,etp2=stp1,stp2
c = 1

#relocate p1
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    return
canvas()
op1,c = visualize_relocation(etp1,[80000,2000],5,200,c)

#relocate p2
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    draw_pair(op1)
    plt.text(90000-1500,27000, 'X', size=30)
    return
op2,c = visualize_relocation(etp2,[100000,2000],5,200,c)

#recombine
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    #draw_pair(op1)
    plt.text(90000-1500,27000, 'X', size=30)
    #draw_pair(op2)
    return

pair1=op1
pair2=op2
rl1 = sorted(r.sample(range(int(pair1[0][0][2][0])+900,int(pair1[0][-1][2][1])-900),4))
rl2 = sorted(r.sample(range(int(pair2[0][0][2][0])+900,int(pair2[0][-1][2][1])-900),4))
for i,i2 in zip(rl1,rl2):
    c = visualize_crossover(pair1,pair2,i,i2,30,c)
    pair1 = crossover(pair1,i)
    pair2 = crossover(pair2,i2)
    c+=1
    
#segregate
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    draw_pair(pair1)
    plt.text(90000-1500,27000, 'X', size=30)
    draw_pair(pair2)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    return

f1,c = visualize_segregation(pair1,pair2,[28500,31500],42000,0.2,200,c)

#relocate f1_1
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    return
canvas()
f1_b1,c = visualize_relocation(f1,[80000,2000],5,200,c)

#relocate f1_2
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    draw_pair(f1_b1)
    plt.text(90000-1500,27000, 'X', size=30)
    return
f1_b2,c = visualize_relocation(f1,[100000,2000],5,200,c)


#recombine
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    #draw_pair(f1_b1)
    plt.text(90000-1500,27000, 'X', size=30)
    return
    return

pair1=f1_b1
pair2=f1_b2
rl1 = sorted(r.sample(range(int(pair1[0][0][2][0])+900,int(pair1[0][-1][2][1])-900),10))
rl2 = sorted(r.sample(range(int(pair2[0][0][2][0])+900,int(pair2[0][-1][2][1])-900),10))
for i,i2 in zip(rl1,rl2):
    c = visualize_crossover(pair1,pair2,i,i2,30,c)
    pair1 = crossover(pair1,i)
    pair2 = crossover(pair2,i2)
    c+=1
    
#segregate into f2
def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    draw_pair(pair1)
    plt.text(90000-1500,27000, 'X', size=30)
    draw_pair(pair2)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    plt.plot(30000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    return

f2,c = visualize_segregation(pair1,pair2,[28500,31500],26500,0.2,200,c)

#show 7 f2 progeny
f2_1 = make_progeny(f1,f1,50)
f2_2 = make_progeny(f1,f1,50)
f2_3 = make_progeny(f1,f1,50)
f2_5 = make_progeny(f1,f1,50)
f2_6 = make_progeny(f1,f1,50)
f2_7 = make_progeny(f1,f1,50)

def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    #draw_pair(pair1)
    #plt.text(90000-1500,27000, 'X', size=30)
    #draw_pair(pair2)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    plt.plot(30000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f2)
    plt.subplot().add_patch(patches.Rectangle((27000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(30000-3100, 26500-3000, "F2 Plant 4")
    
    draw_pair_at_location(f2_1,4500-1500,4500+1500,26500,1)
    plt.plot(4500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((4500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(4500-3100, 26500-3000, "F2 Plant 1")
    
    draw_pair_at_location(f2_2,13000-1500,13000+1500,26500,1)
    plt.plot(13000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((13000-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(13000-3100, 26500-3000, "F2 Plant 2")
    
    draw_pair_at_location(f2_3,21500-1500,21500+1500,26500,1)
    plt.plot(21500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((21500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(21500-3100, 26500-3000, "F2 Plant 3")
    
    draw_pair_at_location(f2_5,38500-1500,38500+1500,26500,1)
    plt.plot(38500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((38500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(38500-3100, 26500-3000, "F2 Plant 5")
    
    draw_pair_at_location(f2_6,47000-1500,47000+1500,26500,1)
    plt.plot(47000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((47000-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(47000-3100, 26500-3000, "F2 Plant 6")
    
    plt.text(53500-3100, 26500-1500, "...", size=15)
    
    draw_pair_at_location(f2_7,55500-1500,55500+1500,26500,1)
    plt.plot(55500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((55500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(55000-3100, 26500-3000,"F2 Plant n")
    
    return

canvas()
c+=1
plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
     
#Zoom in on 7 f2 progeny
#show effect of advcancing generation
fn_1 = scale_pair(f2_1,5) 
fn_2 = scale_pair(f2_2,5) 
fn_3 = scale_pair(f2_3,5) 
fn_4 = scale_pair(f2,5) 
fn_5 = scale_pair(f2_5,5) 
fn_6 = scale_pair(f2_6,5) 
fn_7 = scale_pair(f2_7,5) 

for i in range(7):
    fn_1 = make_progeny(fn_1,fn_1,10)
    fn_2 = make_progeny(fn_2,fn_2,10)
    fn_3 = make_progeny(fn_3,fn_3,10)
    fn_4 = make_progeny(fn_4,fn_4,10)
    fn_5 = make_progeny(fn_5,fn_5,10)
    fn_6 = make_progeny(fn_6,fn_6,10)
    fn_7 = make_progeny(fn_7,fn_7,10)
    
    canvas()
    
    draw_pair_at_location(fn_1,60000+4500-1500,60000+4500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+4500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+4500-2100, 57000-3000,"Plant 1")
    
    draw_pair_at_location(fn_2,60000+13000-1500,60000+13000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+13000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+13000-2100, 57000-3000,"Plant 2")
    
    draw_pair_at_location(fn_3,60000+21500-1500,60000+21500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+21500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+21500-2100, 57000-3000,"Plant 3")
    
    draw_pair_at_location(fn_4,60000+30000-1500,60000+30000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+30000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+30000-2100, 57000-3000,"Plant 4")
    
    draw_pair_at_location(fn_5,60000+38500-1500,60000+38500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+38500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+38500-2100, 57000-3000,"Plant 5")
    
    draw_pair_at_location(fn_6,60000+47000-1500,60000+47000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+47000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+47000-2100, 57000-3000,"Plant 6")
    
    plt.text(60000+ 53500-3100, 2000-1500, "...", size=15)
    
    draw_pair_at_location(fn_7,60000+55500-1500,60000+55500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+55500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+55500-2100, 57000-3000,"Plant n")
    
    plt.text(62500,62000, 'Repeated self-crossing removes all heterozygosity.'.format(i+1), size = 16)
    plt.text(80000,58000, 'Generation: F{0}'.format(i+2), size = 20)
    
    c+=1
    plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
    

#show fn progeny on the left
canvas()
    
draw_pair_at_location(fn_1,60000+4500-1500,60000+4500+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+4500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+4500-2100, 57000-3000,"Plant 1")

draw_pair_at_location(fn_2,60000+13000-1500,60000+13000+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+13000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+13000-2100, 57000-3000,"Plant 2")

draw_pair_at_location(fn_3,60000+21500-1500,60000+21500+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+21500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+21500-2100, 57000-3000,"Plant 3")

draw_pair_at_location(fn_4,60000+30000-1500,60000+30000+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+30000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+30000-2100, 57000-3000,"Plant 4")

draw_pair_at_location(fn_5,60000+38500-1500,60000+38500+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+38500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+38500-2100, 57000-3000,"Plant 5")

draw_pair_at_location(fn_6,60000+47000-1500,60000+47000+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+47000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+47000-2100, 57000-3000,"Plant 6")

plt.text(60000+ 53500-3100, 2000-1500, "...", size=15)

draw_pair_at_location(fn_7,60000+55500-1500,60000+55500+1500,2000,1)
plt.subplot().add_patch(patches.Rectangle((60000+55500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(60000+55500-2100, 57000-3000,"Plant n")

plt.text(62500,62000, 'Repeated self-crossing removes all heterozygosity.'.format(i+1), size = 16)
plt.text(80000,58000, 'Generation: F{0}'.format(i+2), size = 20)


plt.plot(30000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot( 4500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(13000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(21500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(38500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(47000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(55500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(30000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot( 4500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(13000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(21500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(38500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(47000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(55500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(30000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot( 4500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(13000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(21500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(38500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(47000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
plt.plot(55500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')


draw_pair_at_location(fn_1,4500-1500,4500+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((4500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(4500-3100, 4000-3000, "F8 Plant 1")

draw_pair_at_location(fn_2,13000-1500,13000+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((13000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(13000-3100, 4000-3000, "F8 Plant 2")

draw_pair_at_location(fn_3,21500-1500,21500+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((21500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(21500-3100, 4000-3000, "F8 Plant 3")

draw_pair_at_location(fn_4,30000-1500,30000+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((27000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(30000-3100, 4000-3000, "F8 Plant 4")

draw_pair_at_location(fn_5,38500-1500,38500+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((38500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(38500-3100, 4000-3000, "F8 Plant 5")

draw_pair_at_location(fn_6,47000-1500,47000+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((47000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(47000-3100, 4000-3000, "F8 Plant 6")

plt.text(53500-3100, 4000-1500, "...", size=15)

draw_pair_at_location(fn_7,55500-1500,55500+1500,4000,0.2)
plt.subplot().add_patch(patches.Rectangle((55500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
plt.text(55000-3100, 4000-3000,"F8 Plant n")

c+=1
plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)

    
#clear right side    
f8_1 = fn_1
f8_2 = fn_2
f8_3 = fn_3
f8_4 = fn_4
f8_5 = fn_5
f8_6 = fn_6
f8_7 = fn_7

def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    #plt.grid()
    plt.subplot().get_xaxis().set_visible(False)
    plt.subplot().get_yaxis().set_visible(False)
    plt.title("RILs and NILs: Repeat Inbred Lines and Near Isogenic Lines", size=20)
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    draw_pair(etp1)
    draw_pair(etp2)
    plt.text(30000-1000,61500, 'X', size=20)
    #draw_pair(pair1)
    #plt.text(90000-1500,27000, 'X', size=30)
    #draw_pair(pair2)
    plt.plot(30000-400, 55000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f1)
    plt.plot(30000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    draw_pair(f2)
    plt.subplot().add_patch(patches.Rectangle((27000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(30000-3100, 26500-3000, "F2 Plant 4")
    
    draw_pair_at_location(f2_1,4500-1500,4500+1500,26500,1)
    plt.plot(4500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((4500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(4500-3100, 26500-3000, "F2 Plant 1")
    
    draw_pair_at_location(f2_2,13000-1500,13000+1500,26500,1)
    plt.plot(13000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((13000-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(13000-3100, 26500-3000, "F2 Plant 2")
    
    draw_pair_at_location(f2_3,21500-1500,21500+1500,26500,1)
    plt.plot(21500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((21500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(21500-3100, 26500-3000, "F2 Plant 3")
    
    draw_pair_at_location(f2_5,38500-1500,38500+1500,26500,1)
    plt.plot(38500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((38500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(38500-3100, 26500-3000, "F2 Plant 5")
    
    draw_pair_at_location(f2_6,47000-1500,47000+1500,26500,1)
    plt.plot(47000-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((47000-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(47000-3100, 26500-3000, "F2 Plant 6")
    
    plt.text(53500-3100, 26500-1500, "...", size=15)
    
    draw_pair_at_location(f2_7,55500-1500,55500+1500,26500,1)
    plt.plot(55500-400, 40000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.subplot().add_patch(patches.Rectangle((55500-3000, 25000), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(55000-3100, 26500-3000,"F2 Plant n")
    
    
    plt.plot(30000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot( 4500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(13000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(21500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(38500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(47000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(55500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(30000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot( 4500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(13000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(21500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(38500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(47000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(55500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(30000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot( 4500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(13000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(21500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(38500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(47000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
    plt.plot(55500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')


    draw_pair_at_location(f8_1,4500-1500,4500+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((4500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(4500-3100, 4000-3000, "F8 Plant 1")

    draw_pair_at_location(f8_2,13000-1500,13000+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((13000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(13000-3100, 4000-3000, "F8 Plant 2")

    draw_pair_at_location(f8_3,21500-1500,21500+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((21500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(21500-3100, 4000-3000, "F8 Plant 3")

    draw_pair_at_location(f8_4,30000-1500,30000+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((27000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(30000-3100, 4000-3000, "F8 Plant 4")

    draw_pair_at_location(f8_5,38500-1500,38500+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((38500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(38500-3100, 4000-3000, "F8 Plant 5")

    draw_pair_at_location(f8_6,47000-1500,47000+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((47000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(47000-3100, 4000-3000, "F8 Plant 6")

    plt.text(53500-3100, 4000-1500, "...", size=15)

    draw_pair_at_location(f8_7,55500-1500,55500+1500,4000,0.2)
    plt.subplot().add_patch(patches.Rectangle((55500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(55000-3100, 4000-3000,"F8 Plant n")

    return


canvas()
c+=1
plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)

    
    
    
#NILs

#define 7 loci on 7 plants
#largest red loci from the top that isn't at the ends
fn_1 = save_pair_at_location(fn_1,fn_1[0][0][1][0],fn_1[1][0][1][0],2000,1)
fn_2 = save_pair_at_location(fn_2,fn_2[0][0][1][0],fn_2[1][0][1][0],2000,1)
fn_3 = save_pair_at_location(fn_3,fn_3[0][0][1][0],fn_3[1][0][1][0],2000,1)
fn_4 = save_pair_at_location(fn_4,fn_4[0][0][1][0],fn_4[1][0][1][0],2000,1)
fn_5 = save_pair_at_location(fn_5,fn_5[0][0][1][0],fn_5[1][0][1][0],2000,1)
fn_6 = save_pair_at_location(fn_6,fn_6[0][0][1][0],fn_6[1][0][1][0],2000,1)
fn_7 = save_pair_at_location(fn_7,fn_7[0][0][1][0],fn_7[1][0][1][0],2000,1)

l1 = largest_nontelomeric_loci_from_top(fn_1,'red')
l2 = largest_nontelomeric_loci_from_top(fn_2,'red')
l3 = largest_nontelomeric_loci_from_top(fn_3,'red')
l4 = largest_nontelomeric_loci_from_top(fn_4,'red')
l5 = largest_nontelomeric_loci_from_top(fn_5,'red')
l6 = largest_nontelomeric_loci_from_top(fn_6,'red')
l7 = largest_nontelomeric_loci_from_top(fn_7,'red')

qtl1e1 = call_qtl(l1,500,'red')
qtl2e1 = call_qtl(l2,500,'red')
qtl3e1 = call_qtl(l3,500,'red')
qtl4e1 = call_qtl(l4,500,'red')
qtl5e1 = call_qtl(l5,500,'red')
qtl6e1 = call_qtl(l6,500,'red')
qtl7e1 = call_qtl(l7,500,'red')


for i in range(20):
    
    #backcross and selfcross the same plants until defined loci is fixed. 
    fn_1 = select_progeny_from_backcross(fn_1,op2,qtl1e1,'red',500,10)
    fn_2 = select_progeny_from_backcross(fn_2,op2,qtl2e1,'red',500,10)
    fn_3 = select_progeny_from_backcross(fn_3,op2,qtl3e1,'red',500,10)
    fn_4 = select_progeny_from_backcross(fn_4,op2,qtl4e1,'red',500,10)
    fn_5 = select_progeny_from_backcross(fn_5,op2,qtl5e1,'red',500,10)
    fn_6 = select_progeny_from_backcross(fn_6,op2,qtl6e1,'red',500,10)
    fn_7 = select_progeny_from_backcross(fn_7,op2,qtl7e1,'red',500,10)
    
    for i2 in range(8):
        fn_1 = make_progeny(fn_1,fn_1,10)
        fn_2 = make_progeny(fn_2,fn_2,10)
        fn_3 = make_progeny(fn_3,fn_3,10)
        fn_4 = make_progeny(fn_4,fn_4,10)
        fn_5 = make_progeny(fn_5,fn_5,10)
        fn_6 = make_progeny(fn_6,fn_6,10)
        fn_7 = make_progeny(fn_7,fn_7,10)
    
    print("Made it through a cycle.")
    
    canvas()
    
    draw_pair_at_location(fn_1,60000+4500-1500,60000+4500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+4500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+4500-2100, 57000-3000,"Plant 1")
    
    draw_pair_at_location(fn_2,60000+13000-1500,60000+13000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+13000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+13000-2100, 57000-3000,"Plant 2")
    
    draw_pair_at_location(fn_3,60000+21500-1500,60000+21500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+21500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+21500-2100, 57000-3000,"Plant 3")
    
    draw_pair_at_location(fn_4,60000+30000-1500,60000+30000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+30000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+30000-2100, 57000-3000,"Plant 4")
    
    draw_pair_at_location(fn_5,60000+38500-1500,60000+38500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+38500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+38500-2100, 57000-3000,"Plant 5")
    
    draw_pair_at_location(fn_6,60000+47000-1500,60000+47000+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+47000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+47000-2100, 57000-3000,"Plant 6")
    
    plt.text(60000+ 53500-3100, 2000-1500, "...", size=15)
    
    draw_pair_at_location(fn_7,60000+55500-1500,60000+55500+1500,2000,1)
    plt.subplot().add_patch(patches.Rectangle((60000+55500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
    plt.text(60000+55500-2100, 57000-3000,"Plant n")
    
    plt.text(62500,62000, 'Repeated selection of backcross progeny allows fine mapping of QTLs.'.format(i+1), size = 12)
    plt.text(72000,58000, 'Generation: Purebred F8_BC{0}'.format(i+1), size = 20)
    
    c+=1
    plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
    

# #show fn progeny on the left
# canvas()
    
# draw_pair_at_location(fn_1,60000+4500-1500,60000+4500+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+4500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+4500-2100, 57000-3000,"Plant 1")

# draw_pair_at_location(fn_2,60000+13000-1500,60000+13000+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+13000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+13000-2100, 57000-3000,"Plant 2")

# draw_pair_at_location(fn_3,60000+21500-1500,60000+21500+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+21500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+21500-2100, 57000-3000,"Plant 3")

# draw_pair_at_location(fn_4,60000+30000-1500,60000+30000+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+30000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+30000-2100, 57000-3000,"Plant 4")

# draw_pair_at_location(fn_5,60000+38500-1500,60000+38500+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+38500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+38500-2100, 57000-3000,"Plant 5")

# draw_pair_at_location(fn_6,60000+47000-1500,60000+47000+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+47000-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+47000-2100, 57000-3000,"Plant 6")

# plt.text(60000+ 53500-3100, 2000-1500, "...", size=15)

# draw_pair_at_location(fn_7,60000+55500-1500,60000+55500+1500,2000,1)
# plt.subplot().add_patch(patches.Rectangle((60000+55500-3000, 500), 6000, 53000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(60000+55500-2100, 57000-3000,"Plant n")

# plt.text(62500,62000, 'Repeated self-crossing removes all heterozygosity.'.format(i+1), size = 16)
# plt.text(80000,58000, 'Generation: F{0}'.format(i+2), size = 20)


# plt.plot(30000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot( 4500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(13000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(21500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(38500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(47000-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(55500-400, 20000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(30000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot( 4500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(13000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(21500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(38500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(47000-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(55500-400, 19000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(30000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot( 4500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(13000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(21500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(38500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(47000-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')
# plt.plot(55500-400, 18000, markersize=20,marker=r'$\downarrow$', color = 'black')


# draw_pair_at_location(fn_1,4500-1500,4500+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((4500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(4500-3100, 4000-3000, "F8 Plant 1")

# draw_pair_at_location(fn_2,13000-1500,13000+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((13000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(13000-3100, 4000-3000, "F8 Plant 2")

# draw_pair_at_location(fn_3,21500-1500,21500+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((21500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(21500-3100, 4000-3000, "F8 Plant 3")

# draw_pair_at_location(fn_4,30000-1500,30000+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((27000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(30000-3100, 4000-3000, "F8 Plant 4")

# draw_pair_at_location(fn_5,38500-1500,38500+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((38500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(38500-3100, 4000-3000, "F8 Plant 5")

# draw_pair_at_location(fn_6,47000-1500,47000+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((47000-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(47000-3100, 4000-3000, "F8 Plant 6")

# plt.text(53500-3100, 4000-1500, "...", size=15)

# draw_pair_at_location(fn_7,55500-1500,55500+1500,4000,0.2)
# plt.subplot().add_patch(patches.Rectangle((55500-3000, 4000-1500), 6000, 13000, linewidth=1, edgecolor='black', facecolor='none'))
# plt.text(55000-3100, 4000-3000,"F8 Plant n")

# c+=1
# plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
