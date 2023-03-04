
'''
    OLR_tracker.py is designed to identify and track regions highlighted in
    binary maps. It was originally designed to track the Madden-Julian
    Oscillation (MJO) in NOAA Interpolated OLR but can theoretically be used
    to track anything as long as two conditions are met.
        1. The binary maps fed into the tracking algorithm highlight the
           phenomena/object one wants to track
        2. The speed of the object relative to its size on the binary maps
           is small enough such that on subsequent frames

    The main tracking function blob_track takes 3 inputs:
        1. bin_map this is array of binary maps you want to track the designated
           object in. The tracking algorithm assumes a specific orientiation,
           namely (time or Z, lat (y), lon(x)) and for best results binary maps
           should be of this form
        2. min_blob_size the minimum size you expect a tracked object to be on
           the binary maps you input as an integer. Must be greater than or
           equal to 1. Setting this parameter to a larger value will improve the
           speed of the algorithm.
        3. rank_match, set to False by default, is a boolean variable. By default
           tries to match potential continuations of tracked objects based on
           distance (e.g. An animal that temporarily runs behind a tree then
           reappears on the other side). It is also possible to match objects
           based on the direction of travel. If there are multiple possible
           continuations of a tracked object instead of choosing the closest
           possible match rank_match will rank them by distance and direction
           of propagation. The ranking is equal between the two categories.
           If two objects tie in the ranking system, e.g an object is ranked 1st
           for distance and 2nd for angle and another is ranked 1st for angle
           and 2nd distance similarities in propagation speed (velocity) will
           be used to break the tie.

    The output of the tracking algorithm is a Python Dictionary of tracked
    objects. The algorithm is agnostic of tracked object properties other than
    size and when comparing possible continuations to increase its
    generalizability. Because of this the output dictionary should be filtered
    using where tracked objects start, end, propagate, how fast they move, the
    direction in which they move, size, angle of propagation and anything that
    can be calculated from those quantities. With this one should be able to
    refine the results to the exact objects one wishes to identify.
'''


'''
    IMPORTS GO HERE
'''
import numpy as np


'''
    CLASS GOES HERE
'''

class Blob:
    def __init__(self,points,init_frame):
        self.points = sorted(list(set(points)))
        xps = []
        yps = []
        for i in range(len(points)):
            xps.append(points[i][0])
            yps.append(points[i][1])
        self.xpoints = xps
        self.ypoints = yps
        self.cx = 0
        self.cy = 0
        self.vx = 0
        self.vy = 0
        #store old values here
        self.old_points = [self.points]
        self.old_xp = [self.xpoints]
        self.old_yp = [self.ypoints]
        self.old_vx = [0]
        self.old_vy = [0]
        self.old_cx = []
        self.old_cy = []
        self.centroid_calc()
        self.dead = False
        #an array containing the index number of frames the blob existed on
        self.frames = [init_frame]

    def centroid_calc(self):
        x_vals = self.xpoints
        y_vals = self.ypoints

        x_cent = np.mean(x_vals)
        y_cent = np.mean(y_vals)
        self.cx = x_cent
        self.cy = y_cent
        self.old_cx.append(x_cent)
        self.old_cy.append(y_cent)
        if len(self.old_cx) > 1:
            self.vx = self.old_cx[-1] - self.old_cx[-2]
            self.vy = self.old_cy[-1] - self.old_cy[-2]
            self.old_vx.append(self.vx)
            self.old_vy.append(self.vy)



    def update(self,new_ps,frame_num):
        self.points = sorted(list(set(new_ps)))
        xps = []
        yps = []
        for i in range(len(new_ps)):
            xps.append(new_ps[i][0])
            yps.append(new_ps[i][1])

        self.xpoints = xps
        self.ypoints = yps
        self.old_xp.append(self.xpoints)
        self.old_yp.append(self.ypoints)
        self.old_points.append(self.points)
        self.frames.append(frame_num)
        self.centroid_calc()


    def declare_dead(self):
        self.dead = True

    def get_mean_angle(self):
        #returns an array of angles for every day of propagation where one is
        #possible to calculate (stationary days are not included) using
        #the circular mean
        #en.wikipedia.org/Circular_mean
        angle_ray = []
        for i in range(len(self.old_vx)):
            if self.old_vx[i] == 0 and self.old_vy[i] > 0:
                angle_ray.append(90)
            elif self.old_vx[i] == 0 and self.old_vy[i] < 0:
                angle_ray.append(270)
            elif self.old_vx[i] > 0 and self.old_vy[i] == 0:
                angle_ray.append(0)
            elif self.old_vx[i] < 0 and self.old_vy[i] == 0:
                angle_ray.append(180)
            elif self.old_vx[i] == 0 and self.old_vy[i] == 0:
                str = 'do nothing'
            else:
                ang = np.abs(np.rad2deg(np.arctan(self.old_vy[i]/self.old_vx[i])))
                if self.old_vx[i] < 0 and self.old_vy[i] > 0:
                    angle_ray.append(ang + 90)
                elif self.old_vx[i] < 0 and self.old_vy[i] < 0:
                    angle_ray.append(ang + 180)
                elif self.old_vx[i] > 0 and self.old_vy[i] < 0:
                    angle_ray.append(ang + 270)
                else:
                    angle_ray.append(ang)
        if len(angle_ray) > 0:
            inv_len = 1./len(angle_ray)
            #now calculate the circular mean
            mean_angle = np.arctan2(inv_len * np.sum(np.sin(np.deg2rad(np.array(angle_ray)))), inv_len * np.sum(np.cos(np.deg2rad(np.array(angle_ray)))))
            return np.rad2deg(mean_angle)
        else:
            return 9999 #nothing should match with this

    def frame_merge(self):
        #if a blob has been combined with another further on in the algorithm
        #I need to handle duplicate frames. Basically using my methodology
        #an object can be in 2 places at once. I can either take the mean of the
        #2 locations to represent the 'true' location, possibly saying the
        #center is absent of the thing you're tracking or use the connection as
        #the new location losing the information about the final frames of the
        #blob that was connected to the new one. I choose option 2.

        #need to make a copy of each stat, a new list for each and only put the
        #ones I want into the new list
        nop = [] #new old points
        noxp = [] #new old x-points
        noyp = [] #new old y-points
        novx = [] #new old velocity, x-direction
        novy = [] #new old velocity, y-direction
        nocx = [] #new old centroid, x-position
        nocy = [] #new old centroid, y-position
        nf = [] #new frames

        #loop through flames to find duplicates that appear AFTER the original
        for i in range(len(self.frames)):
            dup_count = 0
            if i != (len(self.frames) - 1):
                for j in range(i+1,len(self.frames)):
                    if self.frames[i] == self.frames[j]:
                        dup_count += 1
            if dup_count == 0:
                nop.append(self.old_points[i])
                noxp.append(self.old_xp[i])
                noyp.append(self.old_yp[i])
                novx.append(self.old_vx[i])
                novy.append(self.old_vy[i])
                nocx.append(self.old_cx[i])
                nocy.append(self.old_cy[i])
                nf.append(self.frames[i])


        self.old_points = [ele[1] for ele in sorted(zip(nf,nop))]
        self.old_xp = [ele[1] for ele in sorted(zip(nf,noxp))]
        self.old_yp = [ele[1] for ele in sorted(zip(nf,novy))]
        self.old_vx = [ele[1] for ele in sorted(zip(nf,novx))]
        self.old_vy = [ele[1] for ele in sorted(zip(nf,novy))]
        self.old_cx = [ele[1] for ele in sorted(zip(nf,nocx))]
        self.old_cy = [ele[1] for ele in sorted(zip(nf,nocy))]
        self.frames = sorted(nf)



'''
    FUNCTIONS GO HERE
'''

def connect_8(spx,spy,array):
    '''
    inputs:
        spx: Starting Pixel X position, int, the x value for the pixel I am searching around
        spy: Starting Pixel Y position, int, same as spx but for Y
        array: The array I am searching, array

    '''

    #given the starting x and y of a pixel search for connected pixels to define the blob
    con_x = []
    con_y = []
    #I think I gotta have these 8 conditionals, yuck
    if spx+1 < array.shape[0]:
        if array[spx+1,spy] == 1:
            con_x.append(spx+1)
            con_y.append(spy)
    if spx+1 < array.shape[0] and spy+1 < array.shape[1]:
        if array[spx+1,spy+1] == 1:
            con_x.append(spx+1)
            con_y.append(spy+1)
    if spx+1 < array.shape[0] and spy-1 > 0:
        if array[spx+1,spy-1] == 1:
            con_x.append(spx+1)
            con_y.append(spy-1)
    if spy+1 < array.shape[1]:
        if array[spx,spy+1] == 1:
            con_x.append(spx)
            con_y.append(spy+1)
    if spy-1 > 0:
        if array[spx,spy-1] == 1:
            con_x.append(spx)
            con_y.append(spy-1)
    if spx - 1 > 0 and spy+1 < array.shape[1]:
        if array[spx-1,spy+1] == 1:
            con_x.append(spx-1)
            con_y.append(spy+1)
    if spx - 1 > 0:
        if array[spx-1,spy] == 1:
            con_x.append(spx-1)
            con_y.append(spy)
    if spx - 1 > 0 and spy - 1 > 0:
        if array[spx-1,spy-1] == 1:
            con_x.append(spx-1)
            con_y.append(spy)

    return con_x,con_y


def blob_identify(spx,spy,frame):
    '''
        uses 8-connectivity to find the rest of the binary map blob
        spx,spy = starting coordinates
        frame = the image frame I am searching for blobs
    '''
    spx_con,spy_con = connect_8(spx,spy,frame) #get all points connected to the initial guess
    points_checked = [(spx,spy)] #points checked
    points_found = [(spx,spy)] #points found

    for i in range(len(spx_con)): #add all the connected points to the found array
        points_found.append((spx_con[i],spy_con[i]))

    #now I want to check the points that haven't been checked for 8-connectivity
    while(True):
        points_found = sorted(list(set(points_found))) #remove duplicates
        points_checked = sorted(list(set(points_checked)))
        if points_found == points_checked:
            break
        else:
            for i in range(len(points_found)):
                if points_found[i] not in set(points_checked):
                    p_x,p_y = connect_8(points_found[i][0],points_found[i][1], frame)
                    points_checked.append((points_found[i][0],points_found[i][1]))
                    for i in range(len(p_x)):
                        points_found.append((p_x[i],p_y[i]))




    return sorted(list(set(points_found)))


def blob_search(frame, st):
    '''
        Given a frame goes through and finds all the blobs in it, if they're above
        a certain size, adds them to a blob dictionary

        st = size threshold, minimum number of points to be considered a blob
    '''
    f_swap = np.swapaxes(frame,0,1) #now frame[x,y] is what I'd expect
    #the the square root of the st and round it down
    stsq = int(np.sqrt(st))
    #now make a square just slightly bigger than this to search
    #this almost minimizes rows/columns searched and guarantees that the gaps in between searched rows/columns
    #are smaller than the ST so nothing will be missed
    #a rectangular versus square shape variation is technically better
    #because you can theoretically vary that by +- 1 to get just below the
    #ideal blobject size something to implement in a future version once I think
    #of a clever way to do it
    arr_x = np.arange(0,f_swap.shape[0],stsq+1)
    arr_y = np.arange(0,f_swap.shape[1],stsq+1)
    pf = []
    blob_dict = {}


    for i in range(len(arr_x)):
        for j in range(len(arr_y)): #these two together loop through the array
            tup_x = arr_x[i]
            tup_y = arr_y[j]
            tup_ch = (tup_x,tup_y) #check if the point at tup_ch has already been found
            if (tup_ch not in pf) and (f_swap[arr_x[i],arr_y[j]] == 1):
                #if the point hasn't already been found
                #use 8-connectivity on the point to get everything it is attached to
                bps = blob_identify(arr_x[i],arr_y[j],f_swap) #points of the blob
                for k in range(len(bps)): #add all the found points to pf
                    pf.append(bps[k])
                pf = sorted(list(set(pf))) #remove any duplicates in pf
                if len(bps) >= st: #check if the blob is above a certain threshhold
                    blob_dict[len(blob_dict)] = bps



    return blob_dict


def bst_calc(bm_point,bd_point):
    # gets a blob_similarity value for two blobs to determine if they're the same object

    bin_ray = np.zeros(len(bd_point))
    for i in range(len(bd_point)):
        if bd_point[i] in set(bm_point):
            bin_ray[i] = 1

    bs_val = np.sum(bin_ray) / len(bin_ray)

    return bs_val


def blob_merge_check(blob_dict,blob_master,blob_similarity_threshold,frame_num):
    '''
        this checks if there have been any mergers between blobs, it checks how many
        blobs from previous frame have a sufficient BST with blobs in a new frame
        if two blobs meet BST threshold the older one will become dominant and the
        smaller blob will be declared dead

        blob_dict is the dictionary of blob points found with blob search
        blob_master is the dictionary of previously identified blobjects
        blob_similarity_threshold, what percentage of the new points must overlap with
            previously identified points to declare the two blobs are the same
    '''

    #so I want to loop through the current frames blob_dict and check how many active blobs match
    #up with a blob_dict blob
    for i in range(len(blob_dict)):
        bte = blob_dict[i] #list of points for the blob to examine
        sim_blobs = 0
        sim_blob_inds = []
        #now I want to loop through the blob_master blobs
        for j in range(len(blob_master)):
            if blob_master[j].dead == False:
                blob_sim = bst_calc(blob_master[j].points,bte)
                if blob_sim > blob_similarity_threshold:
                    sim_blobs += 1
                    sim_blob_inds.append(j)

        #the cases are, sim blobs = 0 so its just a new blob
        # sim blobs = 1, its a previously known blob
        # sim blobs > 1, two blobs merged, I only care about this case
        # If the blob isn't bigger I need to declare it dead since it has been
        # absorbed by the larger blob
        if sim_blobs > 1:
            older_blob_age = -9999
            older_blob_ind = -9999
            for K in range(len(sim_blob_inds)):
                if len(blob_master[sim_blob_inds[K]].frames) > older_blob_age:
                    older_blob_age = len(blob_master[sim_blob_inds[K]].frames)
                    older_blob_ind = K
                else:
                    blob_master[sim_blob_inds[K]].declare_dead()
    



def blob_checker(blob_dict,blob_master, blob_similarity_threshold,frame_num):
    '''
        this checks if any blobs have died or new blobs have been added,
        if yes it will update the blobs

        blob_dict is the dictionary of blob points found with blob search
        blob_master is the dictionary of previously identified blobjects
        blob_similarity_threshold, what percentage of the new points must overlap with
            previously identified points to declare the two blobs are the same
    '''
    #first thing to do is do my check for merged blobs
    blob_merge_check(blob_dict,blob_master,blob_similarity_threshold,frame_num)

    #binary array to mark blob_dict blobs that have been checked 0 = not checked, 1 = checked
    blob_dict_bin = np.zeros(len(blob_dict))
    mb_count = 0
    #loop through the blob_master dict and look for blobs that aren't marked as dead
    for i in range(len(blob_master)):
        if blob_master[i].dead == False:
            #so the blob isn't marked as dead, now let's see if there are any matching blobs
            bs_test_ray = np.zeros(len(blob_dict_bin))
            for j in range(len(blob_dict_bin)):
                if blob_dict_bin[j] == 0:
                    bs_test = bst_calc(blob_master[i].points,sorted(list(set(blob_dict[j]))))
                    if bs_test >= blob_similarity_threshold:
                        bs_test_ray[j] = 1
            # so there are now three outcomes, 1 new blob matches an old blob
            # multiple new blobs match an old blob, the new blob doesn't match any
            # the 1:1 case
            if np.sum(bs_test_ray) == 1:
                bi = np.where(bs_test_ray == 1)[0][0]
                blob_master[i].update(sorted(list(set(blob_dict[bi]))),frame_num)
                blob_dict_bin[bi] = 1
            elif np.sum(bs_test_ray) > 1: #multiple blobs have a decent bst
                #only take the oldest blob
                mb_count += 1
                bis = np.where(bs_test_ray == 1)[0]
                max_len_ind = -9999
                max_len = -9999
                for k in range(len(bis)):
                    lc = len(blob_dict[bis[k]])
                    blob_dict_bin[bis[k]] = 1
                    if lc > max_len:
                        max_len = lc
                        max_len_ind = bis[k]
                blob_master[i].update(sorted(list(set(blob_dict[max_len_ind]))),frame_num)
            else: #the case where no blobs match the existing blob
                blob_master[i].declare_dead()

    #now that we've checked all the existing blobs add blobs that didn't match any exisitng blobs
    # to blob master
    for i in range(len(blob_dict_bin)):
        if blob_dict_bin[i] == 0:
            blob_master[len(blob_master)] = Blob(sorted(list(set(blob_dict[i]))),frame_num)

    return blob_master, mb_count

def distance_form(x1,y1,x2,y2):
    #distance formula, nothing fancy
    d = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    return d

def angle_dist_calc(blob1,blob2):
    alpha1 = blob1.get_mean_angle()
    alpha2 = blob2.get_mean_angle()
    #doing it this way guarantees the smaller angular distance
    abs_dif = np.abs(alpha1 - alpha2)
    abs_dif_min = np.abs(alpha1 - alpha2 - 360)

    return np.min([abs_dif,abs_dif_min])

def blob_uniter(blob_master, rank_match = False):
    '''
        Look through the blob_master array to find blobs that are the same one just broken up
        by the inherent variability of the MJO over the MC region
        - Look at the ending location and frame of each blob
        - see if any other blob that initiates within 10 frames and 30 degree radius lasts at least 10 days
        - if so connect it as a continuation of the ended blob
        -give users the option of preferring a ranked angle + distance match
         where blobs are ranked by disance and angle similarity and the lowest
         combined rank blob is used
    '''
    united_blob_master = {} #updated blob master array that has united all blobs
    united_inds = [] #indexes of blobs that were united with another blob
    for i in range(len(blob_master)):
        ef = blob_master[i].frames[-1] #ending frame
        ecx = blob_master[i].cx #ending x-centroid position
        ecy = blob_master[i].cy #ending y-centroid position
        update_bool = False
        matching_inds = []
        for J in range(len(blob_master)):
            if J != i:
                #check if initial frame is within 10 of ending frame and the blob lasts 3 days
                if blob_master[J].frames[0] - ef < 10 and blob_master[J].frames[0] - ef > -5:
                    if len(blob_master[J].frames) >= 2: #also check to make sure it isn't too far N or S
                        if np.abs(distance_form(ecx,ecy,blob_master[J].old_cx[0],blob_master[J].old_cy[0])) <= 15 and np.abs(ecy - blob_master[J].old_cy[0]) <= 3:
                            if (ef > 15 and blob_master[J].frames[0] > 5) or (ef < 15):
                                #now I can assume they're the same blob
                                #I want to add the points of this blob to the original blob
                                matching_inds.append(J) #see how many blobs match
        if len(matching_inds) == 1:
            J = matching_inds[0]
            for K in range(len(blob_master[J].frames)):
                blob_master[i].update(blob_master[J].old_points[K],blob_master[J].frames[K])
            united_inds.append(J)
            united_blob_master[len(united_blob_master)] = blob_master[i]
            update_bool = True
        elif len(matching_inds) > 1: #multiple blobs could be united with the one I'm tracking
            #let's take only the closest
            cent_dist = 9999
            cent_ind = 9999
            if rank_match == False:
                for M in range(len(matching_inds)):
                    if np.abs(distance_form(ecx,ecy,blob_master[matching_inds[M]].old_cx[0],blob_master[matching_inds[M]].old_cy[0])) < cent_dist:
                        cent_dist = np.abs(distance_form(ecx,ecy,blob_master[matching_inds[M]].old_cx[0],blob_master[matching_inds[M]].old_cy[0]))
                        cent_ind = M
                J = matching_inds[cent_ind]
            elif rank_match == True:
                rank_sums = np.zeros(len(matching_inds))
                dists = []
                angle_dif = []
                #get the mean angle of the blob I am trying to connect
                for M in range(len(matching_inds)):
                    dists.append(np.abs(distance_form(ecx,ecy,blob_master[matching_inds[M]].old_cx[0],blob_master[matching_inds[M]].old_cy[0])))
                    angle_dif.append(np.abs(angle_dist_calc(blob_master[i],blob_master[matching_inds[M]])))

                #zip and sort
                #make a copy of the inds of the matching inds
                zip_copy_mi_angle = list(np.arange(0,len(matching_inds),1))
                zip_copy_mi_dist = list(np.arange(0,len(matching_inds),1))
                zcmi_angle_ranked = sorted(zip(angle_dif,zip_copy_mi_angle))
                zcmi_dist_ranked = sorted(zip(dists,zip_copy_mi_dist))
                #I want the rank of each matching blob from 0 -> N
                #this gives me the indexes listed in order of rank
                angle_ranks = [zcmi_angle_ranked[ele][1] for ele in range(len(zcmi_angle_ranked))]
                dist_ranks = [zcmi_dist_ranked[ele][1] for ele in range(len(zcmi_dist_ranked))]
                for rank in range(len(angle_ranks)):
                    rank_sums[angle_ranks[rank]] += rank
                    rank_sums[dist_ranks[rank]] += rank

                #now get the minimum rank sum
                min_rank_sum = np.nanmin(rank_sums)
                rank_sub = rank_sums - min_rank_sum #this is an array where 0 is the best option
                if len(np.where(rank_sub == 0)) < 2:
                    J = np.where(rank_sub == 0)[0][0]
                elif len(np.where(rank_sub == 0)) >= 2:
                    comp_blob_ave_v = np.sign(np.mean(blob_master[i].old_vx)) * np.sqrt(np.mean(blob_master[i].old_vx)**2 + np.mean(blob_master[i].old_vy)**2)
                    v_dif = 9999
                    con_blob_ind = 9999
                    for blob in range(len(np.where(rank_sub == 0))):
                        dif_calc = np.abs(comp_blob_ave_v - np.sign(np.mean(blob_master[np.where(rank_sub == 0)[blob]].old_vx**2 + np.mean(blob_master[np.where(rank_sub == 0)[blob]].old_vy**2))))
                        if dif_calc < v_dif:
                            v_dif = dif_calc
                            con_blob_ind = blob
                    J = np.where(rank_sub == 0)[con_blob_ind][0]

            for K in range(len(list(blob_master[J].frames))):
                blob_master[i].update(blob_master[J].old_points[K],blob_master[J].frames[K])
            united_inds.append(J)
            united_blob_master[len(united_blob_master)] = blob_master[i]
            united_blob_master[len(united_blob_master)-1].frame_merge()
            update_bool = True

        if update_bool == False: #no continuation blob was found so just put what I originally had in here
            if i in set(united_inds): #looking at a blob that was already united
                action_str = 'do nothing'
            else:
                united_blob_master[len(united_blob_master)] = blob_master[i]


    return united_blob_master

def blob_track(bin_map, min_blob_size):
    '''
        Given a binary map of form (time,lat,lon) it tracks
        blobs using the blob-tracking algorithm above
    '''
    blob_master = {}
    for i in range(bin_map.shape[0]):
        if i == 0: #initialize everything
            #get a dictionary of found blobs
            blob_dict = blob_search(bin_map[i],min_blob_size)
            #now use this to populate blob_master
            for J in range(len(blob_dict)):
                blob_master[len(blob_master)] = Blob(blob_dict[J],0)
        else:
            #use the functions I use to update blob_master
            blob_dict = blob_search(bin_map[i], min_blob_size)
            blob_master,mbc = blob_checker(blob_dict,blob_master,.1,i)
    
    #run blob_uniter over the whole thing
    blob_master = blob_uniter(blob_master, rank_match = True)
    #return the master blob dictionary
    for i in range(len(blob_master)):
        blob_master[i].frame_merge()
    return blob_master
