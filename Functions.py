import random as r
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe


#segment = (color,start,end)
#segments = [segment,segment] = [(color,start,end),(color,start,end)]
#pair = [segments,segments] = [[segment,segment],[segment,segment]] = [[(color,start,end),(color,start,end)],[(color,start,end),(color,start,end)]] 

def draw_segment(color, start, end):
    plt.plot(start, end, marker = '', linewidth = 6, path_effects=[mpe.Stroke(linewidth=7, foreground='black'),mpe.Normal()],
        color = color)
    return

def draw_segments(segments):
    for color,start,end in segments:
        draw_segment(color,start,end)
    return

def draw_pair(pair):
    sgs1,sgs2 = pair
    draw_segments(sgs1)
    draw_segments(sgs2)
    return

def draw_pair_at_location(pair,x1,x2,y,scale):
    sgs1,sgs2 = pair
    yo=y
    for color,start,end in sgs1:
        draw_segment(color,[x1,x1],[y,y+((end[1]-end[0])*scale)])
        y = y+((end[1]-end[0])*scale)
    y=yo
    for color,start,end in sgs2:
        draw_segment(color,[x2,x2],[y,y+((end[1]-end[0])*scale)])
        y = y+((end[1]-end[0])*scale)
        
def save_pair_at_location(pair,x1,x2,y,scale):
    sgs1,sgs2 = pair
    nsgs1,nsgs2 = [],[]
    yo=y
    for color,start,end in sgs1:
        nsgs1.append((color,[x1,x1],[y,y+((end[1]-end[0])*scale)]))
        y = y+((end[1]-end[0])*scale)
    y=yo
    for color,start,end in sgs2:
        nsgs2.append((color,[x2,x2],[y,y+((end[1]-end[0])*scale)]))
        y = y+((end[1]-end[0])*scale)
    return [nsgs1,nsgs2]
        
def scale_pair(pair,scale):
    sgs1,sgs2 = pair
    nsgs1,nsgs2 = [],[]
    y = sgs1[0][2][0]
    for color,start,end in sgs1:
        nsgs1.append((color,start,[y,y+((end[1]-end[0])*scale)]))
        y = y+((end[1]-end[0])*scale)
    y = sgs1[0][2][0]
    for color,start,end in sgs2:
        nsgs2.append((color,start,[y,y+((end[1]-end[0])*scale)]))
        y = y+((end[1]-end[0])*scale)
    return [nsgs1,nsgs2]



def sgs_to_df(sgs):
    df = pd.DataFrame(sgs)
    df['s1'] = df[1].apply(lambda x: x[0]).astype(int)
    df['s2'] = df[1].apply(lambda x: x[1]).astype(int)
    df['e1'] = df[2].apply(lambda x: x[0]).astype(int)
    df['e2'] = df[2].apply(lambda x: x[1]).astype(int)
    return df[[0,'s1','s2','e1','e2']]
def df_to_sgs(df):
    sgs = []
    for i,r in df.iterrows():
        sgs.append((r[0],[r['s1'],r['s2']],[r['e1'],r['e2']]))
    return sgs
def mergebreaks(df):
    t = []
    i=0
    #while len(df) > df[0].nunique():
    y = df[0]
    while sum(y.groupby((y != y.shift()).cumsum()).cumcount()) != 0: #while there are consecutive colors
        #print(sum(y.groupby((y != y.shift()).cumsum()).cumcount()),i)
        if df.loc[i,0] == df.iloc[i+1,0]:
            df.loc[i,'e2'] = df.loc[i+1,'e2']
            df = df.drop(i+1).reset_index(drop=1)
            #print(len(df))
            y = df[0]
        else:
            i+=1
    return df


#take a pair, execute a crossover at a location and return the pair
def crossover(pair,location):
    sgs1,sgs2=pair
    #location = location+sgs1[0][2][0]
    nsgs1,nsgs2=[],[]
    
    df1 = sgs_to_df(sgs1)
    df2 = sgs_to_df(sgs2)
    s1,s2=sgs1[0][1],sgs2[0][1]
    
    sgs1p1 = df1[df1['e1'] < location]
    sgs1p1.iloc[-1,-1] = location
    sgs1p2 = df1[df1['e2']> location]
    sgs1p2.iloc[0,3] = location
    
    sgs2p1 = df2[df2['e1'] < location]
    sgs2p1.iloc[-1,-1] = location
    sgs2p2 = df2[df2['e2']> location]
    sgs2p2.iloc[0,3] = location
    
    td1 = pd.concat([sgs1p1,sgs2p2],axis = 0).reset_index(drop=1)
    td2 = pd.concat([sgs2p1,sgs1p2],axis = 0).reset_index(drop=1)
    td1['s1'],td1['s2'] = s1
    td2['s1'],td2['s2'] = s2
    
    td1,td2 = mergebreaks(td1),mergebreaks(td2)

    return [df_to_sgs(td1),df_to_sgs(td2)]

#takes a pair of pairs, execute n crossovers in n locations and return a pair of pairs
def meiosis(pair1,pair2,co_count):
    rl1 = sorted(r.sample(range(int(pair1[0][0][2][0])+900,int(pair1[0][-1][2][1])-900),co_count))
    rl2 = sorted(r.sample(range(int(pair2[0][0][2][0])+900,int(pair2[0][-1][2][1])-900),co_count))
    for i in rl1:
        pair1 = crossover(pair1,i)
    for i in rl2:
        pair2 = crossover(pair2,i)
    return [pair1, pair2]

#take two pairs and return a pair
def segregate(pair1,pair2):
    c1,c2= r.sample([0,1],1)[0], r.sample([0,1],1)[0]
    return [pair1[c1],pair2[c2]]

def make_progeny(pair1,pair2,co_count):
    t1,t2 = meiosis(pair1,pair2,co_count)
    return segregate(t1,t2)


def change_starts(pair,new_start_left,new_start_right):
    sgs1,sgs2 = pair
    nsgs1,nsgs2 = [],[]
    for i in range(len(sgs1)):
        nsgs1.append((sgs1[i][0], [new_start_left,new_start_left], sgs1[i][2]))
    for i in range(len(sgs2)):
        nsgs2.append((sgs2[i][0], [new_start_right,new_start_right], sgs2[i][2]))
    return [nsgs1,nsgs2]

def change_y_scale_up(pair,new_y,upscale):
    sgs1,sgs2 = pair
    nsgs1,nsgs2 = [],[]
    last_increase = 0
    for i in range(len(sgs1)):
        #change_y
        ns1 = (sgs1[i][0], sgs1[i][1], [new_y+(sgs1[i][2][0]-sgs1[0][2][0]),new_y+(sgs1[i][2][0]-sgs1[0][2][0])+(sgs1[i][2][1]-sgs1[i][2][0])])
        #rescale
        y_start = ns1[2][0] + last_increase
        new_length = (ns1[2][1]-ns1[2][0])*upscale
        last_increase = last_increase+new_length-(ns1[2][1]-ns1[2][0])
        ns1 = (ns1[0],ns1[1],[y_start,y_start+new_length])
        nsgs1.append(ns1)
    last_increase = 0
    for i in range(len(sgs2)):
        #change_y
        ns1 = (sgs2[i][0], sgs2[i][1], [new_y+(sgs2[i][2][0]-sgs2[0][2][0]),new_y+(sgs2[i][2][0]-sgs2[0][2][0])+(sgs2[i][2][1]-sgs2[i][2][0])])
        #rescale
        y_start = ns1[2][0] + last_increase
        new_length = (ns1[2][1]-ns1[2][0])*upscale
        last_increase = last_increase+new_length-(ns1[2][1]-ns1[2][0])
        ns1 = (ns1[0],ns1[1],[y_start,y_start+new_length])
        nsgs2.append(ns1)
        
    return [nsgs1,nsgs2]

def visualize_relocation(pair,destination,upscale,steps,c):
    sgs1,sgs2 = pair
    center = (sgs1[0][1][0] + sgs2[0][1][0])/2
    
    #c=0
    for step in range(steps+1):
        current_x = center + step*((destination[0]-center)/steps)
        current_y = sgs1[0][2][0] + step*((destination[1]-sgs1[0][2][0])/steps)
        
        left = current_x-(center-sgs1[0][1][0])
        right = current_x+(sgs2[0][1][0]-center)
        
        scale = 1+step*((upscale-1)/steps)
        nsgs1,nsgs2 = change_y_scale_up(pair,current_y,scale)
        
        nsgs1,nsgs2 = change_starts([nsgs1,nsgs2],left,right)
        
        canvas()
        draw_pair([nsgs1,nsgs2])
        c+=1
        plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
        
    return [nsgs1,nsgs2],c

def visualize_segregation(pair1,pair2,destination_x,destination_y,downscale,steps,c):
    c1,c2= r.sample([0,1],1)[0], r.sample([0,1],1)[0]
    interchromosomal_distance = 3000
    sgs1,sgs2 = pair1[c1],pair2[c2]
    pair = [sgs1,sgs2]
    left_x = sgs1[0][1][0]
    right_x = sgs2[0][1][0]
    y = sgs1[0][2][0]
    destination_x_left,destination_x_right = destination_x
    
    #c=0
    for step in range(steps+1):
        current_left_x = left_x + step*((destination_x_left-left_x)/steps)
        current_right_x = right_x + step*((destination_x_right-right_x)/steps)
        current_y = y + step*((destination_y-y)/steps)
        
        scale = 1+step*((downscale-1)/steps)
        nsgs1,nsgs2 = change_y_scale_up(pair,current_y,scale)
        
        nsgs1,nsgs2 = change_starts([nsgs1,nsgs2],current_left_x,current_right_x)
        
        canvas()
        draw_segments(nsgs1)
        draw_segments(nsgs2)
        c+=1
        plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
    
    return [nsgs1,nsgs2],c

#only adjusts upper segment start when crossover is at junction
#issue. when upper segment length < 900, the start of the next segment protrudes downwards. #fixed. I think.

def visualize_crossover(pair1,pair2,location1,location2,steps,c):
    adj=900
    adj2=450
    flag=False
    sgs1,sgs2 = pair1
    left,right = sgs1[0][1][0],sgs2[0][1][0]
    #c=0
    for step in range(steps+1):
        if (location1 <= sgs1[0][2][0]+adj) or (location1 >= sgs1[-1][2][1]-adj):
            return #crossing over within adj range at the ends results in truncation.
        
        canvas()
        sgs1,sgs2 = pair1
        location = location1
        #draw sgs1   
        x_left = sgs1[0][1][0]
        for i in range(len(sgs1)):
            if (sgs1[i][2][1] > location) and (sgs1[i][2][0] < location): #are we in the segment that needs to be broken?
                #break it into two.
                lower_segment=(sgs1[i][0],sgs1[i][1],[sgs1[i][2][0],location])
                #draw lower pair
                draw_segment(lower_segment[0],lower_segment[1], [lower_segment[2][0], lower_segment[2][1]-0])
                plt.plot([sgs1[0][1][0],sgs2[0][1][0]],[location+adj2,location+adj2], linestyle = '--', color = 'green')
                #update x coordinates and draw upper pair
                x_left = sgs1[i][1][0] + step*((right-left)/steps)
                upper_segment=(sgs1[i][0],[x_left, x_left],[location,sgs1[i][2][1]])
                #when upper segment length < adj, open a flag
                if (upper_segment[2][1] - location) < 900:
                    flag = True
                else:
                    draw_segment(upper_segment[0],upper_segment[1],[upper_segment[2][0]+adj, upper_segment[2][1]])
            elif (sgs1[i][2][0] == location):
                x_left = sgs1[i][1][0] + step*((right-left)/steps)
                draw_segment(sgs1[i][0],[x_left,x_left],[sgs1[i][2][0]+adj,sgs1[i][2][1]])
                plt.plot([sgs1[0][1][0],sgs2[0][1][0]],[location+adj2,location+adj2], linestyle = '--', color = 'green')
            else:
                #draw the rest
                if flag:
                    flag = False
                    draw_segment(sgs1[i][0],[x_left,x_left],[location+adj,sgs1[i][2][1]])
                    print('something')
                else:
                    draw_segment(sgs1[i][0],[x_left,x_left],sgs1[i][2])
                
        #draw sgs2
        x_right = sgs2[0][1][0]
        for i in range(len(sgs2)):
            if (sgs2[i][2][1] > location) and (sgs2[i][2][0] < location): #are we in the segment that needs to be broken?
                #break it into two.
                lower_segment=(sgs2[i][0],sgs2[i][1],[sgs2[i][2][0],location])
                #draw lower pair
                draw_segment(lower_segment[0],lower_segment[1],[lower_segment[2][0], lower_segment[2][1]-0])
                #update x coordinates and draw upper pair
                x_right = sgs2[i][1][0] - step*((right-left)/steps)
                upper_segment=(sgs2[i][0],[x_right,x_right],[location,sgs2[i][2][1]])
                #when upper segment length < adj, open a flag
                if (upper_segment[2][1] - location) < 900:
                    flag = True
                else:
                    draw_segment(upper_segment[0],upper_segment[1],[upper_segment[2][0]+adj, upper_segment[2][1]])
            elif (sgs2[i][2][0] == location):
                x_right = sgs2[i][1][0] - step*((right-left)/steps)
                draw_segment(sgs2[i][0],[x_right,x_right],[sgs2[i][2][0]+adj,sgs2[i][2][1]])
            else:
                #draw the rest
                if flag:
                    flag = False
                    draw_segment(sgs2[i][0],[x_right,x_right],[location+adj,sgs2[i][2][1]])
                    print('something')
                else:
                    draw_segment(sgs2[i][0],[x_right,x_right],sgs2[i][2])
        
        sgs1,sgs2 = pair2
        location = location2
        #draw sgs1   
        x_left = sgs1[0][1][0]
        for i in range(len(sgs1)):
            if (sgs1[i][2][1] > location) and (sgs1[i][2][0] < location): #are we in the segment that needs to be broken?
                #break it into two.
                lower_segment=(sgs1[i][0],sgs1[i][1],[sgs1[i][2][0],location])
                #draw lower pair
                draw_segment(lower_segment[0],lower_segment[1], [lower_segment[2][0], lower_segment[2][1]-0])                
                plt.plot([sgs1[0][1][0],sgs2[0][1][0]],[location+adj2,location+adj2], linestyle = '--', color = 'green')
                #update x coordinates and draw upper pair
                x_left = sgs1[i][1][0] + step*((right-left)/steps)
                upper_segment=(sgs1[i][0],[x_left, x_left],[location,sgs1[i][2][1]])
                #when upper segment length < adj, open a flag
                if (upper_segment[2][1] - location) < 900:
                    flag = True
                else:
                    draw_segment(upper_segment[0],upper_segment[1],[upper_segment[2][0]+adj, upper_segment[2][1]])
            elif (sgs1[i][2][0] == location):
                x_left = sgs1[i][1][0] + step*((right-left)/steps)
                draw_segment(sgs1[i][0],[x_left,x_left],[sgs1[i][2][0]+adj,sgs1[i][2][1]])
                plt.plot([sgs1[0][1][0],sgs2[0][1][0]],[location+adj2,location+adj2], linestyle = '--', color = 'green')
            else:
                #draw the rest
                if flag:
                    flag = False
                    draw_segment(sgs1[i][0],[x_left,x_left],[location+adj,sgs1[i][2][1]])
                    print('something')
                else:
                    draw_segment(sgs1[i][0],[x_left,x_left],sgs1[i][2])
                
        #draw sgs2
        x_right = sgs2[0][1][0]
        for i in range(len(sgs2)):
            if (sgs2[i][2][1] > location) and (sgs2[i][2][0] < location): #are we in the segment that needs to be broken?
                #break it into two.
                lower_segment=(sgs2[i][0],sgs2[i][1],[sgs2[i][2][0],location])
                #draw lower pair
                draw_segment(lower_segment[0],lower_segment[1],[lower_segment[2][0], lower_segment[2][1]-0])
                #update x coordinates and draw upper pair
                x_right = sgs2[i][1][0] - step*((right-left)/steps)
                upper_segment=(sgs2[i][0],[x_right,x_right],[location,sgs2[i][2][1]])
                #when upper segment length < adj, open a flag
                if (upper_segment[2][1] - location) < 900:
                    flag = True
                else:
                    draw_segment(upper_segment[0],upper_segment[1],[upper_segment[2][0]+adj, upper_segment[2][1]])
            elif (sgs2[i][2][0] == location):
                x_right = sgs2[i][1][0] - step*((right-left)/steps)
                draw_segment(sgs2[i][0],[x_right,x_right],[sgs2[i][2][0]+adj,sgs2[i][2][1]])
            else:
                #draw the rest
                if flag:
                    flag = False
                    draw_segment(sgs2[i][0],[x_right,x_right],[location+adj,sgs2[i][2][1]])
                    print('something')
                else:
                    draw_segment(sgs2[i][0],[x_right,x_right],sgs2[i][2]) 
        
        c+=1
        plt.savefig('testf2/a{0}.png'.format(str(c).zfill(5)),format = 'png', dpi=200)
    return c


def canvas():
    plt.figure(figsize=(16,9))
    plt.xlim(0,120000)
    plt.ylim(0,68000)
    plt.grid()
    plt.axvline(x = 60000, color = 'black', linestyle = '-')
    plt.text(30000-1000,61500, 'X', size=20)
    return