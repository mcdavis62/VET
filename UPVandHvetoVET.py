"""
Created on Thursday October 6 15:14:00 2022

Last modified on Friday June 16 2023

@author: Michael C. Davis

Function Listing:

    - UPV_scraper:
        query UPV daily result files
        makes the UPV.txt, UPV_merge.txt, UPV_segs.txt, and UPV_segs.txt for all of the data and for each round
        
    - UPV_sort_on_col:
        sorts a data array based on the values in a single column
        makes the UPV_sorted.txt file containing the veto times in ascending order
        
    - query_UPV_flag:
        query the UPV daily result files
        makes the path to the UPV_sorted.txt file
        
    - hveto:
        query hveto daily result files
        makes the hveto.txt, hveto_merge.txt, and hveto_segs.txt files for all of the data and for each round
        
    - query_hveto_flag:
        query hveto daily result files
        makes the path to the hveto_sorted.txt file
        also recreates the paths for the hveto.txt, hveto_merge.txt, and hveto_segs.txt files
        
    - hveto_sort_on_col:
        sorts a data array based on the values in a single column
        makes the UPV_sorted.txt file containing the veto times in ascending order
        
    - make_flag_xml:
        makes condidate flag xml file for all of the data and then the round winners
        
    - query_cwb_events:
        query cwb event file from a CIT cluster
        creates the cwb.lcf and the EVENTS.txt files
        
    - Write_config:
        creates the cwb_daily.ini file
        
    - make_excecutable:
        makes the run.sh file (executable)
        
    - Wrapper
        assembles the files needed to execute daily UPV and Hveto segments against the cWB results


"""



#=================================================================================================================
#=================================================================================================================

''' Start of UPV_scraper() '''

def UPV_scraper(year, month, day, dest_dir, ifo, verbose = False):
    
    '''
    description:
                query UPV daily result files
                makes the UPV.txt, UPV_merge.txt, UPV_segs.txt, and UPV_segs.txt for all of the
                    data and for each round
                    
    :param year: year (int)
    :param month: month (int)
    :param day: day (int)
    :param dest_dir: destination directory (str)
    :param ifo: detector, e.g. H1 or L1 (str)
    
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    if verbose:
        print(f"UPV")
    
    if verbose:
        print(f"Processing date: {year}-{month}-{day}")
    
    # Converts the year, month, day to GPS date
    gps = to_gps(datetime.datetime(year, month, day))
    date = gps.gpsSeconds
    date_end = date + 86400
    
    if verbose:
        print(f"GPS time: {date}-{date_end}\n")

    # Path to the UPV veto segment round winners page
    directory = f"/home/detchar/public_html/NewUPV/O3/results/O3-H-{date}_{date_end}/"

    # Path for the .html page
    file = directory + "index.html"
    f = open(file, 'r')

    # Counter for the number of channel sources; for pages that have more than one channel
    i = 0

    # Empty array to hold the channel source names
    source_channel = np.array([])
    # Empty array to hold the channel segments
    segs = np.array([])
    
    if verbose:
        print(f"|| Getting round winners...")

    # Reads the .html page to find how many channels there are and their names
    for line in f:
        # If the line starts with "<h2>"
        if line[0:4] == '<h2>':
            # If the 13-17 characters have "S[number]" then that line has the channel name
            if line[12:16] == f'"S{i}"':
                # Splits the line by the spaces
                line_array = line.split(' ')
                # Appends the empty array and adds the 5th item in the split array (the channel source)
                source_channel = np.append(source_channel, line_array[5])
                # Adds to the channel source counter
                i += 1

    # The number of channels
    channels = len(source_channel)
    
    if verbose:
        print(f"-- There are {channels} channels")
        print(f"-- The channels are: ")
    
    # Empty array holders
    flag = np.array([])
    flag_file = np.array([])
    N = 1
    UPV_channel_names = np.array([])
    
    # Formats the channel names for the Flag Name and File variables
    for k in range(len(source_channel)):
        ch = source_channel[k]
        root_ch = ch.split(':')[1]
        flag_ch = f'{ifo}:UPV_ROUND{N}_{root_ch}:1'
        flag = np.append(flag, flag_ch)
        flag_ch_file = f'{ifo}-UPV_ROUND{N}_{root_ch}.xml'
        flag_file = np.append(flag_file, flag_ch_file)
        N += 1
        # The names of the channels
        UPV_channel_names = np.append(UPV_channel_names, root_ch)
        if verbose:
            print(f"\t{root_ch}")
      
    if verbose:
        print(f"-> Done (round winners)")
        
    
    # --------------------------------------------------------+
    # Check if the {ifo}-UPV_ALLROUNDS.xml file already exits    
    # --------------------------------------------------------+
    for j in range(len(source_channel)):
        if os.path.exists(f'{source_channel[j]}.xml'):
            # If it already exits, remove it
            os.remove(f'{source_channel[j]}.xml') 
            
    if verbose:
        print(f"|| Getting UPV flag name...")

    # Creates a dictionary to hold the channel names in the Flag Name variable
    keys = range(len(source_channel))
    
    # Sets the first item in the dictionary to have all of the data
    UPV_flag_name = {
        
        0 : f'{ifo}:UPV_ALLROUNDS:1'
        
    }
    
    # Fills in the dictionary with remaining rounds
    for k in keys:
        UPV_flag_name[k+1] = flag[k]
        
    if verbose:
        print(f"-> Done (UPV flag name)")
        
    if verbose:
        print(f"|| Getting UPV flag file...")
     
    # Creates a dictionary to hold the channel names in the Flag File variable
    # Sets the first item in the dictionary to have all of the data
    UPV_flag_file = {
        
        0 : f'{ifo}-UPV_ALLROUNDS.xml'

    }
    
    # Fills in the dictionary with remaining rounds
    for k in keys:
        UPV_flag_file[k+1] = flag_file[k]
    
    if verbose:
        print(f"-> Done (UPV flag file)")
    
    # Creates empty arrays to hold the start and end times
    start_time = np.array([])
    end_time = np.array([])

    # Path to the UPV_segs.txt file (list of the path of UPV files)
    UPV_dest_path_seg_file = os.path.join(dest_dir, 'UPV_segs.txt')
        
    # -------------------------------------+
    # Check if the UPV_segs.txt file already exits    
    # -------------------------------------+
    if os.path.exists(UPV_dest_path_seg_file):
        # If it already exits, remove it
        os.remove(UPV_dest_path_seg_file) 

    # Path to the UPV_merge.txt file (contains all the contents in all the UPV result files; veto segment, start/end times, veto duration)
    UPV_dest_path_merge = os.path.join(dest_dir, 'UPV_merge.txt')

    # -------------------------------------+
    # Check if the UPV_merge.txt file already exits    
    # -------------------------------------+
    if os.path.exists(UPV_dest_path_merge):
        # If it already exits, remove it
        os.remove(UPV_dest_path_merge) 

    # Path to the UPV.txt file (start and end times of all round UPV results)
    dest_path_UPV = os.path.join(dest_dir, 'UPV.txt')


    '''
    # -------------------------------------+
    # Check if the UPV.txt file already exits    
    # -------------------------------------+
    if os.path.exists(dest_path_UPV):
        # If it already exits, remove it
        os.remove(dest_path_UPV) 
    '''    

    if verbose:
        print(f"|| Creating UPV_segs.txt, UPV_merge.txt, and UPV.txt files for each round...")

    # For every source channel
    for k in range(len(source_channel)):
        # The path to the veto file
        veto_file = directory + f"{ifo}:GDS-CALIB_STRAIN/" + source_channel[k] + "/veto_" + source_channel[k] + ".txt"

        #print(veto_file)
        
        # Path to the UPV_segs.txt file (list of the path of UPV files) for the round winners
        UPV_dest_path_seg_file_winner = os.path.join(dest_dir, f'UPV_segs_rd{k + 1}.txt')

        # Path to the UPV_merge.txt file for the round winners (contains all the contents in all the UPV result files; veto segment, start/end times, veto duration)
        UPV_dest_path_merge_winner = os.path.join(dest_dir, f'UPV_merge_rd{k + 1}.txt')
        
        # Path to the UPV.txt file for the round winners (start and end times of all round UPV results)
        dest_path_UPV_winner = os.path.join(dest_dir, f'UPV_rd{k + 1}.txt')

        seg_round = np.array([veto_file])
        
        # Creates the segs.txt file for the round winners
        np.savetxt(UPV_dest_path_seg_file_winner, seg_round, delimiter = '\n', fmt = '%s')

        # Creates the segs.txt file containing the channel paths
        segs = np.append(segs, veto_file)
        np.savetxt(UPV_dest_path_seg_file, segs, delimiter = '\n', fmt = '%s')

        start, end = np.loadtxt(veto_file, unpack=True)

        # The time duration of each veto
        time_step_winner = end - start

        # The number of veto segments
        veto_segments_winner = np.arange(0, len(start), 1)

        # Creates the UPV_merge.txt file for the round winners containing the veto segment, start time, end time, and the veto duration
        np.savetxt(UPV_dest_path_merge_winner, np.transpose([veto_segments_winner, start, end, time_step_winner]), delimiter = '\t')    

        # Creates the UPV.txt file for round winners
        np.savetxt(dest_path_UPV_winner, np.transpose([start, end]), delimiter = ' ')
        
        if verbose:
            if os.path.exists(f'{dest_dir}/UPV_segs_rd{k+1}.txt') and os.path.exists(f'{dest_dir}/UPV_merge_rd{k+1}.txt') and os.path.exists(f'{dest_dir}/UPV_rd{k+1}.txt'):
                print(f"-- Created UPV_segs.txt, UPV_merge.txt, and UPV.txt files for round {k+1} - {UPV_channel_names[k]}")
            else:
                raise IOError(f"UPV_segs.txt, UPV_merge.txt, and UPV.txt files were NOT created for round {k+1} - {UPV_channel_names[k]}")

        # Appends the start and end time arrays
        start_time = np.append(start_time, start)
        end_time = np.append(end_time, end)

    if verbose:
        print(f"-> Done (UPV_segs.txt, UPV_merge.txt, and UPV.txt files for all rounds)")

    if verbose:
        print(f"|| Creating UPV.txt file for all data...")
        
    # Creates the UPV.txt file containing the start and end times in two columns - for the entire data
    np.savetxt(dest_path_UPV, np.transpose([start_time, end_time]), delimiter = ' ')  
    
    if verbose:
        if os.path.exists(f'{dest_dir}/UPV.txt'):
            print(f"-> Done (UPV.txt)")
        else:
            raise IOError(f"UPV.txt was NOT created")

    # The time duration of each veto
    time_step = end_time - start_time

    # The number of veto segments
    veto_segments = np.arange(0, len(start_time), 1)
    
    if verbose:
        print(f"|| Creating UPV_merge.txt file for all data...")

    # Creates the UPV_merge.txt file containing the veto segment, start time, end time, and the veto duration - for the entire data
    np.savetxt(UPV_dest_path_merge, np.transpose([veto_segments, start_time, end_time, time_step]), delimiter = '\t')

    if verbose:
        if os.path.exists(f'{dest_dir}/UPV_merge.txt'):
            print(f"-> Done (UPV_merge.txt)")
        else:
            raise IOError(f"UPV_merge.txt was NOT created")

    # -------------------------------------------------------+
    # Check if the UPV_seg.txt has zero UPV results    |
    # -------------------------------------------------------+
    if os.stat(UPV_dest_path_seg_file).st_size == 0:
        # print('File is empty')
        raise IOError("no UPV results on this day")

    # Returns the number of channels, the channel names, and the dictionaries for the Flag Name and File    
    return channels, UPV_channel_names, UPV_flag_name, UPV_flag_file

''' End of UPV_scraper() '''

#=================================================================================================================
#=================================================================================================================

''' Start of UPV_sort_on_col() '''

# Creates the UPV_sorted txt file (veto times in ascending order)
def UPV_sort_on_col(data, col = 1, ascending = True, verbose = False):
        
    '''
    description:
                sorts a data array based on the values in a single column
                makes the UPV_sorted.txt file containing the veto times in ascending order
    
    usage: sort_on_col(data, col=0, ascending=True)
    
    :param data: data to sort on columns (array) - data = np.loadtxt(f'{dest_dir}/UPV.txt', unpack = True)
    :param col: column to sort on (int)
    :param ascending: if True, return ascending order; if False, return descending order (bool) - default = True
    
    output:                  
    sorted Data (array)
    
    requires numpy library
    '''
    
    import numpy as np

    
    if not ascending:
        data = data[flipud(data[:,col].argsort())]
    else:
        data = data[data[:,col].argsort()]
        
    return data

''' End of UPV_sort_on_col() '''

#=================================================================================================================
#=================================================================================================================

''' Start of query_UPV_flag() '''

def query_UPV_flag(year, month, day, dest_dir, ifo, verbose = False):
    
    '''
    description:
                query the UPV daily result files.
                makes the path to the UPV_sorted.txt file

    param year: year (int)
    param month: month (int)
    param day: day (int)
    param dest_dir: destination directory (str)
    param ifo: detector, e.g. H1 or L1 (str)
    
    return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    # Converts the year, month, day to GPS date
    gps = to_gps(datetime.datetime(year, month, day))
    date = gps.gpsSeconds
    date_end = date + 86400
    

    # The directory to the UPV results; the date must be in GPS time
    base_dir = f"/home/detchar/public_html/NewUPV/O3/results/O3-H-{date}_{date_end}/{ifo}:GDS-CALIB_STRAIN" # .format(DATE=date)
    
    
    # Path to the UPV_sorted.txt file (sorted start and end times in ascending order)
    dest_path_sorted_UPV = os.path.join(dest_dir, 'UPV_sorted.txt') 
    
    # Path to the UPV_sorted.txt files for all rounds
    
    
    return None

''' End of query_UPV_flag() '''

#=================================================================================================================
#=================================================================================================================

''' Start of hveto() '''

def hveto(year, month, day, dest_dir, ifo, verbose = False):
    
    '''
    description:
                query hveto daily result files
                makes the hveto.txt, hveto_merge.txt, and hveto_segs.txt files for all of the
                    data and for each round

    :param year: year (int)
    :param month: month (int)
    :param day: day (int)
    :param dest_dir: destination directory (str)
    :param ifo: detector, e.g., H1 or L1 (str)
    :return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    if verbose:
        print(f"\nHveto")
    
    if verbose:
        print(f"Processing date: {year}-{month}-{day}")
    
    # Converts the year, month, day to GPS date
    gps = to_gps(datetime.datetime(year, month, day))
    gps = gps.gpsSeconds
    gps_end = gps + 86400
    
    if verbose:
        print(f"GPS time: {gps}-{gps_end}\n")
    
    # The hveto directory requires the format of the date to be a string in the following example: 20190801
    year = str(year)
    
    month = f"{month:02}"
    
    day = f"{day:02}"
    
    date = year + month + day
    
    # Path to the hveto veto segment round winners page
    base_dir = f"/home/detchar/public_html/hveto/day/{date}/latest/segments/"

    # Creates empty arrays to hold the segment, start, end, and time step times
    segment = np.array([])
    start_time = np.array([])
    end_time = np.array([])
    time_step = np.array([])

    # Creates an empty array for the segs variable; later used for the hveto_segs.txt
    segs = np.array([])
    
    if verbose:
        print(f"|| Getting round winners...")

    # Displays the number of veto segment txt files are in the diectory; used for the loop
    cmd = f'ls {base_dir}{ifo}-HVETO_VETO_SEGS_ROUND_*.txt -l | wc -l'
    rounds = os.popen(cmd).read()
    rounds = rounds[0:-1]
    rounds = int(rounds)
    
    # Loads in the summary statistics txt file
    summary = np.loadtxt(f'/home/detchar/public_html/hveto/day/{date}/latest/summary-stats.txt', unpack = True, dtype = str)

    # From the summary statistics, gather the channel names of the round winners
    channels = summary[1]
    
    if verbose:
        print(f"-- There are {len(channels)} channels")
        print(f"-- The channels are: ")
    
    # Empty array holders
    flag = np.array([])
    flag_file = np.array([])
    N = 1
    hveto_channel_names = np.array([])
    
    # Formats the channel names for the Flag Name and File variables
    for k in range(len(channels)):
        ch = channels[k]
        root_ch = ch.split(':')[1]
        flag_ch = f'{ifo}:HVETO_ROUND{N}_{root_ch}:1'
        flag = np.append(flag, flag_ch)
        flag_ch_file = f'{ifo}-HVETO_ROUND{N}_{root_ch}.xml'
        flag_file = np.append(flag_file, flag_ch_file)
        N += 1
        # The names of the channels
        hveto_channel_names = np.append(hveto_channel_names, root_ch)
        if verbose:
            print(f"\t{root_ch}")
            
    if verbose:
        print(f"-> Done (round winners)")
        

    # --------------------------------------------------------+
    # Check if the {ifo}-hveto_ALLROUNDS.xml file already exits    
    # --------------------------------------------------------+
    for j in range(len(channels)):
        if os.path.exists(f'{channels[j]}.xml'):
            # If it already exits, remove it
            os.remove(f'{channels[j]}.xml') 
    
    # The number of round winning channels
    hveto_channels = len(channels)
    
    if verbose:
        print(f"|| Getting hveto flag name...")
    
    # Creates a dictionary to hold the channel names in the Flag Name variable
    keys = range(len(channels))
    
    # Sets the first item in the dictionary to have all of the data
    hveto_flag_name = {
        
        0 : f'{ifo}:HVETO_ALLROUNDS:1'

    }
    
    # Fills the dictionary with the remaining rounds
    for k in keys:
        hveto_flag_name[k+1] = flag[k]
        
    if verbose:
        print(f"-> Done (hveto flag name)")
        
    if verbose:
        print(f"|| Getting hveto flag file...")

    # Creates a dictionary to hold the channel names in the Flag File variable
    # Sets the first item in the dictionary to have all of the data
    hveto_flag_file = {
        
        0 : f'{ifo}-HVETO_ALLROUNDS.xml'
   
    }
    
    # Fills in teh dictionary with the remaining rounds
    for k in keys:
        hveto_flag_file[k+1] = flag_file[k]
        
    if verbose:
        print(f"-> Done (hveto flag file)")
        
        
    ''' Rounds '''


    if verbose:
        print(f"|| Creating hveto_segs.txt, hveto_merge.txt, and hveto.txt files for each round...")    
    
    # For every source channel
    for k in range(1, rounds + 1):
        # The path to the veto file
        file = f'/home/detchar/public_html/hveto/day/{date}/latest/segments/{ifo}-HVETO_VETO_SEGS_ROUND_{k}-{gps}-86400.txt'
        if os.path.exists(file):
            seg, start, end, step = np.loadtxt(file, unpack = True)
        
            # Makes the hveto_merge.txt file for every round winner
            hveto_dest_path_merge_winner = os.path.join(dest_dir, f'hveto_merge_rd{k}.txt')
            np.savetxt(hveto_dest_path_merge_winner, np.transpose([seg, start, end, step]), delimiter = ' ')
        
            # Makes the hveto.txt file for every round winner
            dest_path_hveto_winner = os.path.join(dest_dir, f'hveto_rd{k}.txt')
            np.savetxt(dest_path_hveto_winner, np.transpose([start, end]), delimiter = ' ')
            
            
            # Makes the hveto_segs.txt file for every round winner
            round_segs = np.array([file])
            hveto_dest_path_seg_file_winner = os.path.join(dest_dir, f'hveto_segs_rd{k}.txt')
            np.savetxt(hveto_dest_path_seg_file_winner, round_segs, delimiter = '\n', fmt = '%s')
        
            # Appends the empty arrays
            segment = np.append(segment, seg)
            start_time = np.append(start_time, start)
            end_time = np.append(end_time, end)
            time_step = np.append(time_step, step)
            
            if verbose:
                print(f"-- Created hveto_segs.txt, hveto_merge.txt, and hveto.txt files for round {k} - {hveto_channel_names[k-1]}")
                
    if verbose:
        print(f"-> Done (hveto_segs.txt, hveto_merge.txt, hveto.txt for all rounds)")  
        
        
    ''' All Data '''
        
        
    if verbose:
        print(f"|| Creating hveto_merge.txt file for all data...")

    # Path to the hveto_merge.txt file (contains all the contents in all the hveto result files; veto segment, start/end times, veto duration)
    hveto_dest_path_merge = os.path.join(dest_dir, 'hveto_merge.txt')

    # -------------------------------------+
    # Check if the hveto_merge.txt file already exits    
    # -------------------------------------+
    if os.path.exists(hveto_dest_path_merge):
        # If it already exits, remove it
        os.remove(hveto_dest_path_merge)

    # Creates the hveto_merge.txt file containing the veto segment, start time, end time, and the veto duration
    np.savetxt(hveto_dest_path_merge, np.transpose([segment, start_time, end_time, time_step]), delimiter = '\t')

    if verbose:
        if os.path.exists(hveto_dest_path_merge):
            print(f"-> Done (hveto_merge.txt)")

            
    if verbose:
        print(f"|| Creating hveto.txt file for all data...")    

    # Path to the hveto.txt file (start and end times of all round hveto results)
    dest_path_hveto = os.path.join(dest_dir, 'hveto.txt')

    # Creates the hveto.txt file containing the start and end times in two columns
    np.savetxt(dest_path_hveto, np.transpose([start_time, end_time]), delimiter = ' ')

    if verbose:
        if os.path.exists(dest_path_hveto):
            print(f"-> Done (hveto.txt)")
            
    if verbose:
        print(f"|| Creating hveto_segs.txt file for all data...")


    # Path to the hveto_segs.txt file (list of the path of hveto files)
    hveto_dest_path_seg_file = os.path.join(dest_dir, 'hveto_segs.txt')

    # -------------------------------------+
    # Check if the hveto_segs.txt file already exits    
    # -------------------------------------+
    if os.path.exists(hveto_dest_path_seg_file):
        # If it already exits, remove it
        os.remove(hveto_dest_path_seg_file) 

    # Creates the segs.txt file containing the channel paths
    segs = np.append(segs, file)
    np.savetxt(hveto_dest_path_seg_file, segs, delimiter = '\n', fmt = '%s')
    
    if verbose:
        if os.path.exists(hveto_dest_path_seg_file):
            print(f"-> Done (hveto_segs.txt)")

    # -------------------------------------------------------+
    # Check if the hveto_segs.txt has zero hveto results    |
    # -------------------------------------------------------+
    if os.stat(hveto_dest_path_seg_file).st_size == 0:
        # print('File is empty')
        raise IOError("no hveto results on this day")
        
                
    return hveto_channels, hveto_channel_names, hveto_flag_name, hveto_flag_file

''' End of hveto() '''

#=================================================================================================================
#=================================================================================================================

''' Start of query_hveto_flag() '''

def query_hveto_flag(year, month, day, dest_dir, ifo, verbose = False): 
    '''
    description:
                query hveto daily result files
                makes the path to the hveto_sorted.txt file
                also recreates the paths for the hveto.txt, hveto_merge.txt, and hveto_segs.txt files


    :param date: date, e.g, yyyymmdd
    :param dest_dir: destination directory
    :param ifo: ifo, e.g., H1 or L1
    :return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    # Formats the date to work with the hveto directory
    year = str(year)
    
    month = f"{month:02}"
    
    day = f"{day:02}"
    
    date = year + month + day

    # Path to the veto segments
    base_dir = "/home/detchar/public_html/hveto/day/{DATE}/latest/segments".format(DATE=date)
    
    # Path to the file of sorted and end times
    dest_path_sorted_hveto = os.path.join(dest_dir, 'hveto_sorted.txt')
    
    # Path to the list of the path of hveto files
    hveto_dest_path_seg_file = os.path.join(dest_dir, 'hveto_segs.txt')
    
    # Path to the file that contains all the contents in all the hveto result files
    hveto_dest_path_merge = os.path.join(dest_dir, 'hveto_merge.txt')
    
    # Path to the file of start and end times of all round hveto results
    dest_path_hveto = os.path.join(dest_dir, 'hveto.txt')
    

    # -------------------------------------------------+
    # query the files of Hveto result files           |
    # and to write the path to the files to segs.txt  |
    # -------------------------------------------------+
    cmd = f'ls {base_dir}/{ifo}-HVETO_VETO_SEGS*.txt > {hveto_dest_path_seg_file}'
    
    os.system(cmd)

    
    # ----------------------------------+
    # check if the hveto_merge file exits    
    # ----------------------------------+
    if os.path.exists(hveto_dest_path_merge):
        os.remove(hveto_dest_path_merge)  # if it already exits, remove it

    # -------------------------------------------------------+
    # check if the hveto_seg.txt has at least one hvete result    |
    # -------------------------------------------------------+
    if os.stat(hveto_dest_path_seg_file).st_size == 0:
        #print('File is empty')
        raise IOError("no hveto results on this day")
        
    # ---------------------------------------------------------------+
    # merge all the content in the hveto result files it merge.txt  |
    # ---------------------------------------------------------------+
    else:
        with open(hveto_dest_path_seg_file) as fo:
            for line in fo:
                line = line.rstrip('\n')
                command = f'cat {line} >> {hveto_dest_path_merge}'.format(line=line, hveto_dest_path_merge = hveto_dest_path_merge)
                os.system(command)

        # -----------------------------------------------------------------------+
        # write the start and end times of the actives segments into hveto.txt  |
        # -----------------------------------------------------------------------+
        with open(dest_path_hveto, "w") as f:
            with open(hveto_dest_path_merge) as fi:
                for line in fi:
                    start_time = line.split('\t')[1]
                    end_time = line.split('\t')[2]
                    f.write(start_time)
                    f.write('\t')
                    f.write(end_time)
                    f.write('\n')
                    
    return None

''' End of query_hveto_flag() '''

#=================================================================================================================
#=================================================================================================================

''' Start of hveto_sort_on_col() '''

# Creates the hveto_sorted txt file (veto times in ascending order)
def hveto_sort_on_col(data, col = 1, ascending = True, verbose = False):
        
    '''
    description:
                sorts a data array based on the values in a single column
                makes the UPV_sorted.txt file containing the veto times in ascending order
    
    usage: sort_on_col(data, col=0, ascending=True)
    
    :param data: data to sort on columns (array) - 
    :param col: column to sort on (int)
    :param ascending: if True, return ascending order; if False, return descending order (bool) - default = True
    
    output:                  
    sorted Data (array)
    
    requires numpy library
    '''
    
    import numpy as np
    
    if not ascending:
        data = data[flipud(data[:,col].argsort())]
    else:
        data = data[data[:,col].argsort()]
        
    return data

''' End of hveto_sort_on_col() '''

#=================================================================================================================
#=================================================================================================================

''' Start of make_flag_xml() '''

def make_flag_xml(dest_dir, ifo, year, month, day, flag_name, flag_file, active_file, verbose = False):
    '''
    description:
                makes condidate flag xml file for all of the data and then the round winners
                
    usage: run twice - one for UPV and one for hveto

    :param dest_dir: destination directory (str)
    :param ifo: detector, e.g., H1 or L1 (str)
    :param year: year (int)
    :param month: month (int)
    :param day: day (int)
    :param flag_name: flag name - name of flag name which is expected to be a :1 format (str)
        will be taken from the created dictionaries in either UPV_scraper() or hveto()
    :param flag_file: flag file - name of flag file which is expected to be a .xml format (str)
        will be taken from the created dictionaries in either UPV_scraper() or hveto()
    :param active_file: file which contains start and end times of active segments - UPV_sorted.txt or
        hveto_sorted.txt files (str)
    :return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    if verbose:
        print(f"\t|| Creating flag xml file...")
    

    # Time range of flag definition:
    start = to_gps(datetime.datetime(year, month, day))
    stop = start + 86400

    # File name of active (veto) segments and its delimiter:
    active_delim = '\t'

    
    observing = '{IFO}:DMT-ANALYSIS_READY:1'.format(IFO=ifo)
    try:
        known_segs = DataQualityFlag.query_dqsegdb(observing, start, stop).active
    except RuntimeError:
        print("Please type the following")
        print("$ ligo-proxy-init albert.einstein")
        print("by replacing your user name")
        print("then rerun it")
        sys.exit()

    # === Specify definition times, input files, and flag name below:
    path_active_file = os.path.join(dest_dir, active_file)
    
    if verbose:
        print(f"\t-- Reading active file: {path_active_file}")

    # Name of this flag
    flag_name = '{Flag_name}'.format(Flag_name=flag_name)

    # Name of the segment file:
    xml_file = '{Flag_file}'.format(Flag_file=flag_file)

    # path to the xml file
    path_xml_file = os.path.join(dest_dir, xml_file)
    
    # ----------------------------------+
    # check if the path_xml_file already exits    
    # ----------------------------------+
    if os.path.exists(path_xml_file):
        # if it already exits, remove it
        os.remove(path_xml_file)  
    
    if verbose:
        print(f"\t-- Reading path to xml file: {path_xml_file}")

    # Read the active (veto) segments
    active_segs = np.loadtxt(path_active_file, delimiter= ' ')

    # Get an array for the active_start and active_end of each veto segment
    active_start = [active_segs[i, 0].round(2) for i in range(len(active_segs))]
    active_end = [active_segs[i, 1].round(2) for i in range(len(active_segs))]
    
    if verbose:
        print(f"\t-- Creating data quality flag object...")

    # Create a data quality flag object
    flag = DataQualityFlag(flag_name, active=zip(active_start, active_end), known=known_segs)
    
    if verbose:
        print(f"\t-- Writing flag to xml file")

    # Write flag to file - creates the xml file
    flag.write(path_xml_file, format = 'ligolw')
    
    if verbose:
        print(f"\t-> Done (flag xml file)")
            
    return None

''' End of make_flag_xml() '''

#=================================================================================================================
#=================================================================================================================

''' Start of query_cwb_events() '''

def query_cwb_events(dest_dir, ifo, year, month, day, verbose = False):
    '''
    description:
                query cwb event file from a CIT cluster
                creates the cwb.lcf file
                creates the EVENTS.txt file

    :param dest_dir: destination directory
    :param ifo: ifo, e.g., H1 or L1
    :param year: year
    :param month: month
    :param day: day
    :return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb

    if verbose:
        print(f"|| Gertting cWB events file...")

    # Takes either L1 or H1 without the one - L or H
    ifo_for_cwb = ifo[0]  

    # Converts the year, month, day into GPS date
    start = to_gps(datetime.datetime(year, month, day))
    stop = start + 86400

    # ------------------------------------------------------------------------------------+
    # copy the cwb daily EVENTS.txt file from a CIT cluster to the destination directory |
    # ------------------------------------------------------------------------------------+
    
    # Path to the events file
    # If the gps time is before 1 Oct 2019, in O3a
    if start < 1253923218:
        cwb_base_dir = f"/home/waveburst/public_html/online/O3_LH_BurstLF_ONLINE/POSTPRODUCTION/FOM_daily_{start}-{stop}/plotbin2_cut/data/EVENTS.txt"
    # If the gps time is after 1 Nov 2019, in O3b
    elif start >= 1256601618: 
        cwb_base_dir = f"/home/waveburst/public_html/online/O3b_LH_BurstLF_ONLINE/POSTPRODUCTION/FOM_daily_{start}-{stop}/plotbin2_cut/data/EVENTS.txt"
    else:
        print('There is no data in October (during break in O3).')
        
    if verbose:
        print(f"-- Events.txt file from: {cwb_base_dir}")
    
    cmd = f'scp ldas-pcdev1.ligo.caltech.edu:{cwb_base_dir} {dest_dir}'
    
    os.system(cmd)

    # -------------------------------+
    # this portion creates cwb.lcf  |
    # -------------------------------+
    with open(os.path.join(dest_dir, 'cwb.lcf'), 'w') as f:
        f.write(ifo_for_cwb)
        f.write(' ')
        f.write('CWB')
        f.write(' ')
        f.write(str(start.gpsSeconds))
        f.write(' ')
        f.write(str(stop.gpsSeconds))
        f.write(' ')
        f.write('EVENTS.txt')
        
    # print(cwb_base_dir)
    
    # print(cmd)
    
    if verbose:
        print(f"-> Done (cWB events file)")
    
    return None

''' End of query_cwb_events() '''

#=================================================================================================================
#=================================================================================================================

''' Start of Write_config() '''

def Write_config(dest_dir, config_file, IFO, flag_name1, flag_file1, flag_name2, flag_file2, paddings, UPV_channel_names, hveto_channel_names, verbose = False):
    """
    This function is inspired by the file made by Robert Beda
    (see https://ldas-jobs.ligo-la.caltech.edu/~robert.beda/)
    All the credit goes to him.
    Writes a file to perform a VET run with desired parameters,
    without using DQSEGDB.
    
    description:
                creates the cwb_daily.ini file

    :param dest_dir: destination directory (str)
    :param config_file: the configuraton file, e.g., cwb_daily.ini (str)
    :param IFO: detector, e.g., H1 or L1 (str)
    :param flag_name1: flag name in the UPV flag
    :param flag_file1: name of the UPV flag file, which is expected to be a .xml file
    :param flag_name2: flag name in the hveto flag
    :param flag_file2: name of the hveto flag file, which is expected to be a .xml file
    :param paddings: tupple or list of the paddings, ([pre_pad, post_pad], ....), e.g., ([1, 0], [0, 1], [1, 1])
    :return: None

    usage: Write_config(dest_dir, config_file, IFO, flag_name1, flag_file1, flag_name2, flag_file2, paddings)

    """
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    if verbose:
        print(f"|| Writing config file to perform VET...")

    to_write = []
    # This line is sufficient to tell the summary page software that this is a VET run
    to_write.append('[plugins]\n')
    to_write.append('gwvet.tabs =\n')
    to_write.append('\n')
    
    # Standard link to note a VET issue
    to_write.append('[html]\n')
    to_write.append('issues = https://github.com/gwpy/vet/issues\n')
    to_write.append('\n')
    
    # This defines what flags correspond to different states.
    to_write.append('[states]\n')
    to_write.append('Locked = {IFO}:DMT-DC_READOUT_LOCKED:1\n'.format(IFO=IFO))
    to_write.append('Science = {IFO}:DMT-ANALYSIS_READY:1\n'.format(IFO=IFO))
    to_write.append('\n')
    
    #  This defines the location of the segment database.
    to_write.append('[segment-database]\n')
    to_write.append('url = https://segments.ligo.org\n')
    to_write.append('\n')
    
    # for cWB triggers
    to_write.append('[cwb]\n')
    to_write.append("timecolumn = 'time for {IFO} detector'\n".format(IFO=IFO))
    to_write.append("columns = 'time for {IFO} detector','effective correlated amplitude rho','central frequency','sSNR for {IFO} detector'\n".format(IFO=IFO))
    to_write.append('\n')
    
    # default settings
    to_write.append('[DEFAULT]\n')
    to_write.append('type = veto-flag\n')
    to_write.append('event-channel = {IFO}:GDS-CALIB_STRAIN\n'.format(IFO=IFO))
    to_write.append('event-generator = cwb\n')
    to_write.append('event-format = ascii.cwb\n')
    to_write.append("metrics = 'Deadtime',\n")
    to_write.append("          'Efficiency',\n")
    to_write.append("          'Efficiency/Deadtime',\n")
    to_write.append("          'Efficiency | effective correlated amplitude rho>=6.5',\n")
    to_write.append("          'Efficiency/Deadtime | effective correlated amplitude rho>=6.5',\n")
    to_write.append("          'Efficiency | effective correlated amplitude rho>=7',\n")
    to_write.append("          'Efficiency/Deadtime | effective correlated amplitude rho>=7',\n")
    to_write.append("          'Efficiency | effective correlated amplitude rho>=8',\n")
    to_write.append("          'Efficiency/Deadtime | effective correlated amplitude rho>=8',\n")
    to_write.append("          'Use percentage',\n")
    to_write.append("          'Loudest event by effective correlated amplitude rho'\n")
    to_write.append('\n')
    
    # Creates the tab for UPV - all of the data
    to_write.append('[tab-All_UPV]\n')
    to_write.append(f'parent = UPV\n')
    to_write.append(f'name = All Rounds\n')
    to_write.append(f'flags = {flag_name1[0]}\n')
    to_write.append(f'states = Science\n')
    to_write.append(f'segmentfile = {flag_file1[0]}\n')
    to_write.append('\n')
             
    # Creates the tabs for UPV - each round
    for t in range(len(UPV_channel_names)):
        to_write.append(f'[tab-UPV_Round{t+1}]\n')
        to_write.append(f'parent = UPV\n')
        to_write.append(f'name = Round {t+1} - {UPV_channel_names[t]}\n')
        to_write.append(f'flags = {flag_name1[t+1]}\n')
        to_write.append(f'states = Science\n')
        to_write.append(f'segmentfile = {flag_file1[t+1]}\n')
        to_write.append('\n')
    
    # Creates the tab for Hveto - all of the data
    to_write.append('[tab-All_Hveto]\n')
    to_write.append(f'parent = Hveto\n')
    to_write.append(f'name = All Rounds\n')
    to_write.append(f'flags = {flag_name2[0]}\n')
    to_write.append(f'states = Science\n')
    to_write.append(f'segmentfile = {flag_file2[0]}\n')
    to_write.append(f'\n')

    # Creates the tabs for Hveto - each round
    for t in range(len(hveto_channel_names)):
        to_write.append(f'[tab-Hveto_Round{t+1}]\n')
        to_write.append(f'parent = Hveto\n')
        to_write.append(f'name = Round {t+1} - {hveto_channel_names[t]}\n')
        to_write.append(f'flags = {flag_name2[t+1]}\n')
        to_write.append(f'states = Science\n')
        to_write.append(f'segmentfile = {flag_file2[t+1]}\n')
        to_write.append(f'\n')
    
    if paddings is not None:
        # with padding
        for index_pad in range(len(paddings)):
            pre_pad, post_pad = paddings[index_pad]
            if pre_pad != 0 and post_pad == 0:
                tab_name = 'm{Pre_pad}'.format(Pre_pad=int(pre_pad))
                padding = '(-{Pre_pad},{Post_pad})'.format(Pre_pad=int(pre_pad), Post_pad=int(post_pad))
            elif pre_pad == 0 and post_pad !=0:
                tab_name = 'p{Post_pad}'.format(Post_pad=int(post_pad))
                padding = '({Pre_pad},+{Post_pad})'.format(Pre_pad=int(pre_pad), Post_pad=int(post_pad))
            elif pre_pad !=0 and post_pad !=0:
                tab_name = 'm{Pre_pad}p{Post_pad}'.format(Pre_pad=int(pre_pad), Post_pad=int(post_pad))
                padding = '(-{Pre_pad},+{Post_pad})'.format(Pre_pad=int(pre_pad), Post_pad=int(post_pad))
            pad_name = '-{Pre_pad},+{Post_pad}'.format(Pre_pad=int(pre_pad), Post_pad=int(post_pad))
            to_write.append('[tab-{Tab_name}]\n'.format(Tab_name=tab_name))
            to_write.append('parent = padded\n')
            to_write.append('name = {Pad_name}\n'.format(Pad_name=pad_name))
            to_write.append(f'flags = {flag_name1}\n')
            to_write.append('states = Science\n')
            to_write.append(f'segmentfile = {flag_file1}\n')
            to_write.append('padding = {Padding}\n'.format(Padding=padding))
            to_write.append('\n')
    
    # Joins the appended writing to the config file path - creates cwb_daily.ini
    with open(os.path.join(dest_dir, config_file), 'w') as f:
        f.write("".join(to_write))
        f.close()
        
        
    if verbose:
        print(f"-> Done (config file)")

''' End of Write_config() '''
        
#=================================================================================================================
#=================================================================================================================
        
''' Start of make_excecutable() '''
        
def make_excecutable(dest_dir, exe_file, config_file, year, month, day, verbose = False):
    '''
    description:
                makes the run.sh file (executable)
        
    :param dest_dir: destination directory (str)
    :param exe_file: name of the executable file (str)
    :param config_file: name of the configuration file (str)
    :param year: year (int)
    :param month: month (int)
    :param day: day (int)
    :return: None
    '''
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    if verbose:
        print(f"|| Creating an executable file...")

    # Path to the exe_file (run.sh)
    #path_executable = os.path.join('/home/michael.davis/public_html/O3b/', exe_file)
    path_executable = os.path.join(dest_dir, exe_file)
    
    # Converst the year, month, day to GPS date
    start = to_gps(datetime.datetime(year, month, day))
    stop = start + 86400

    # Fills in the commands that will be in the exe_file (run.sh)
    command = []
    command.append('#!/bin/bash\n')
    command.append('source /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/etc/profile.d/conda.sh\n')
    command.append('conda activate /home/detchar/.conda/envs/ligo-summary-3.9\n')
    command.append('export GWPY_USETEX=0\n')
    command.append('gw_summary gps ')
    command.append('{} '.format(int(start)))
    command.append('{} '.format(int(stop)))
    command.append('--event-cache cwb.lcf -f {config_file} --verbose'.format(config_file=config_file))


    # Joins the appended writing to the exe_file - creates run.sh
    with open(path_executable, 'w') as f:
        f.write(''.join(command))
        
    if verbose:
        print(f"-> Done (executable file)")
        
    return None

''' End of make_excecutable() '''

#=================================================================================================================
#=================================================================================================================

''' Start of Wrapper() '''

def Wrapper(ifo, year, month, day, dest_dir = '.', verbose = False):

    """
    Wrapper() assembles the files needed to execute daily UPV segments against the cWB results

    ifo         = H1 or L1 (str)
    year        = 4-digit year (int)
    month       = 2-digit month (int)
    day         = 2-digit day (int)
    dest_dir    = (optional) destination directory, default = '.'
    verbose     = switch to display progess (boolean), default = False


    ''' assembles the files needed to execute daily UPV and hveto segments against the cWB results '''
    ''' executes VET '''
    
    """
    
    import numpy as np
    import os
    import datetime
    import time
    from gwpy.time import to_gps, from_gps
    from gwpy.segments import DataQualityFlag
    import pdb
    
    """
    check_kinit = input(f"In the terminal, have you entered 'kinit', your albert.einstein password, and 'ligo-proxy-init -k'?\nAnswer yes or no")
    if check_kinit == "yes" or check_kinit == "Yes" or check_kinit == "YES":
        pass
    else:
        print(f"Go to your terminal window and enter 'kinit', your albert.einstein password, and 'ligo-proxy-init -k")
    """
    
    
    ''' Checking the provided dates and parameters are valid '''
    
    if not 1 <= month <= 12:
        raise IOError(f"The month provided is not a valid month")
        
    if not 1 <= day <= 31:
        raise IOError(f"The day provided is not a valid day")
        
    if ifo != 'H1' and ifo != "L1":
        raise IOError(f"The IFO provided is not a valid IFO")

    timer_start = time.time()
    timer = False

    '''
    UPV_flag_name = None
    hveto_flag_name = None

    UPV_flag_file = None
    hveto_flag_file = None
    '''
    config_file = 'cwb_daily.ini'

    exe_file = 'run.sh'

    pads = None
    
    #no_data_counter = 0

    # if we are running over a single day...
    if isinstance(day, int):

        # Converts the year, month, day to GPS date
        gps = to_gps(datetime.datetime(year, month, day))
        date = gps.gpsSeconds
        date_end = date + 86400
        
        if verbose:
            print(f"Processing date: {year}-{month}-{day}\n")
           
        if verbose:
            print(f"GPS time: {date}-{date_end}\n")


        # Path to the UPV veto segment round winners page
        directory = f"/home/detchar/public_html/NewUPV/O3/results/O3-H-{date}_{date_end}/"

        # Path for the html page
        file = directory + "index.html"
        
        # If the path for the UPV data exists, then carry on with the wrapper and executing VET
        if os.path.exists(file):
            #print('path exists!')
            #sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
            # custom modules
            #from auto_vet.utils.utils import query_hevto_flag, make_flag_xml, query_cwb_events, Write_config, make_excecutable # the list of glitch classes in each IFO

            work_dir = os.getcwd()

            if True:#__name__ == '__main__':

                # make paddings
                if pads is not None:
                    paddings = [[int(pad.split(',')[0][1:]), int(pad.split(',')[-1][:-1])] for pad in pads]
                else:
                    paddings = None


                # Converts the date to the correct format
                date = "{}{:02d}{:02d}".format(year, month, day)

                # Checks if the Flag Files already exist for UPV and hveto

                # --------------------------------------------------------+
                # Check if the {ifo}-UPV_ALLROUNDS.xml file already exits    
                # --------------------------------------------------------+
                if os.path.exists(f'{dest_dir}{ifo}-UPV_ALLROUNDS.xml'):
                    # If it already exits, remove it
                    os.remove(f'{dest_dir}{ifo}-UPV_ALLROUNDS.xml') 

                # --------------------------------------------------------+
                # Check if the {ifo}-HVETO_ALLROUNDS.xml file already exits    
                # --------------------------------------------------------+
                if os.path.exists(f'{dest_dir}{ifo}-HVETO_ALLROUNDS.xml'):
                    # If it already exits, remove it
                    os.remove(f'{dest_dir}{ifo}-HVETO_ALLROUNDS.xml') 


                # Make the destination directory if it does not exist yet
                try:
                    os.mkdir(os.path.join(dest_dir))
                except OSError:
                    pass


                ''' UPV '''


                # Runs the UPV Scraper function
                channels, UPV_channel_names, UPV_flag_name, UPV_flag_file = UPV_scraper(year, month, day, dest_dir, ifo)  
                if verbose:
                    print(f"|| Creating the UPV_sorted.txt file for all data...")

                # Runs the UPV_sort_on_col function for the entire data
                ascended_UPV = UPV_sort_on_col(data = np.loadtxt(f'{dest_dir}/UPV.txt', unpack = True), col = 1, ascending = True)
                # Arranges the UPV vetoes in ascending order
                start_time_ascend = ascended_UPV[0]
                end_time_ascend = ascended_UPV[1]
                
                # Creates the UPV_sorted.txt file containing the veto times in ascending order
                np.savetxt(dest_dir + '/UPV_sorted_rd0.txt', np.transpose([start_time_ascend, end_time_ascend]), delimiter = ' ')

                if verbose:
                    if os.path.exists(f'{dest_dir}/UPV_sorted_rd0.txt'):
                        print(f"-> Done (UPV_sorted_rd0.txt)")
                    else:
                        raise IOError("UPV_sorted_rd0.txt was NOT created")

                if verbose:
                    print(f"|| Creating the UPV_sorted.txt file for each round...")

                # Runs the UPV_sort_on_col function for every round
                for i in range(1, channels + 1):
                    ascended_UPV = UPV_sort_on_col(data = np.loadtxt(f'{dest_dir}/UPV_rd{i}.txt', unpack = True), col = 1, ascending = True)
                    # Arranges the UPV vetoes in ascending order
                    start_time_ascend = ascended_UPV[0]
                    end_time_ascend = ascended_UPV[1]
                    # Creates the UPV_sorted.txt file containing the veto times in ascending order
                    np.savetxt(dest_dir + f'/UPV_sorted_rd{i}.txt', np.transpose([start_time_ascend, end_time_ascend]), delimiter = ' ')

                    if verbose:
                        if os.path.exists(f'{dest_dir}/UPV_sorted_rd{i}.txt'):
                            print(f"-- Created UPV_sorted file for round {i} - {UPV_channel_names[i-1]}")
                        else:
                            raise IOError(f"UPV_sorted.txt for round {i} (UPV_sorted_rd{i}.txt) was NOT created")

                        
                ''' Query UPV Flag '''


                # Runs the query_UPV_flag function
                query_UPV_flag(year, month, day, dest_dir, ifo)


                
                ''' Hveto '''


                hveto_channels, hveto_channel_names, hveto_flag_name, hveto_flag_file = hveto(year, month, day, dest_dir, ifo)

                if verbose:
                    print(f"|| Creating the hveto_sorted.txt file for all data...")

                # Runs the hveto_sort_on_col function for the entire data
                ascended_hveto = hveto_sort_on_col(data = np.loadtxt(f'{dest_dir}/hveto.txt', unpack = True), col = 1, ascending = True)
                # Arranges the hveto vetoes in ascending order
                hveto_start_time_ascend = ascended_hveto[0]
                hveto_end_time_ascend = ascended_hveto[1]
                # Creates the hveto_sorted.txt file containing the veto times in ascending order
                np.savetxt(dest_dir + '/hveto_sorted_rd0.txt', np.transpose([hveto_start_time_ascend, hveto_end_time_ascend]), delimiter = ' ')

                if verbose:
                    if os.path.exists(f'{dest_dir}/hveto_sorted_rd0.txt'):
                        print(f"-> Done (hveto_sorted_rd0.txt)")
                    else:
                        raise IOError("hveto_sorted_rd0.txt was NOT created")

                if verbose:
                    print(f"|| Creating the hveto_sorted.txt file for each round...")

                # Runs the hveto_sort_on_col function for every round
                for k in range(1, hveto_channels + 1):
                    ascended_hveto = hveto_sort_on_col(data = np.loadtxt(f'{dest_dir}/hveto_rd{k}.txt', unpack = True), col = 1, ascending = True)
                    # Arranges the hveto vetoes in ascending order
                    hveto_start_time_ascend = ascended_hveto[0]
                    hveto_end_time_ascend = ascended_hveto[1]
                    # Creates the hveto_sorted.txt file containing the veto times in ascending order
                    np.savetxt(dest_dir + f'/hveto_sorted_rd{k}.txt', np.transpose([hveto_start_time_ascend, hveto_end_time_ascend]), delimiter = ' ')

                    if verbose:
                        if os.path.exists(f'{dest_dir}/hveto_sorted_rd{k}.txt'):
                            print(f"-- Created hveto_sorted file for round {k} - {hveto_channel_names[k-1]}")
                        else:
                            raise IOError(f"hveto_sorted.txt for round {k} (hveto_sorted_rd{k}.txt) was NOT created")
                        
                        
                ''' Query Hveto Flag '''


                # Runs the query_hveto_flag function
                query_hveto_flag(year, month, day, dest_dir, ifo)



                ''' Make Flag XML '''

                if verbose:
                    print(f"\n|| Creating flag xml files for UPV...\n")

                # Runs the make_flag_xml function for UPV - for all of the data and all rounds
                for n in UPV_flag_name:
                    make_flag_xml(dest_dir, ifo, year, month, day, flag_name = UPV_flag_name[n], flag_file = UPV_flag_file[n], active_file = f'UPV_sorted_rd{n}.txt')
                    
                    if verbose:
                        if n == 0:
                            if os.path.exists(f'{dest_dir}/{UPV_flag_file[n]}'):
                                print(f"\t-- xml file for all data created - {UPV_flag_file[n]}\n")
                            else:
                                raise IOError(f"xml file for all data was NOT created")
                        else:
                            if os.path.exists(f'{dest_dir}/{UPV_flag_file[n]}'):
                                print(f"\t-- xml file for round {n} created - {UPV_flag_file[n]}\n")
                            else:
                                raise IOError(f"xml file for round {n} ({UPV_flag_file[n]}) was NOT created")

                if verbose:
                    print(f"-> Done (flag xml files for all rounds of UPV)")


                if verbose:
                    print(f"\n|| Creating flag xml files for hveto...\n")

                # Runs the make_flag_xml function for hveto - for all of the data and all rounds
                for m in hveto_flag_name:
                    make_flag_xml(dest_dir, ifo, year, month, day, flag_name = hveto_flag_name[m], flag_file = hveto_flag_file[m], active_file = f'hveto_sorted_rd{m}.txt')
                    if verbose:
                        if m == 0:
                            if os.path.exists(f'{dest_dir}/{hveto_flag_file[m]}'):
                                print(f"\t-- xml file for all data created - {hveto_flag_file[m]}\n")
                            else:
                                raise IOError(f"xml file for all data was NOT created")
                        else:
                            if os.path.exists(f'{dest_dir}/{hveto_flag_file[m]}'):
                                print(f"\t-- xml file for round {m} created - {hveto_flag_file[m]}\n")
                            else:
                                raise IOError(f"xml file for round {m} ({hveto_flag_file[m]}) was NOT created")

                if verbose:
                    print(f"-> Done (flag xml files for all rounds of hveto)")


                ''' Query cWB Events '''

                if verbose:
                    print(f"\n|| Creating cwb.lcf and EVENTS.txt files...")

                # Runs the query_cwb_events function
                query_cwb_events(dest_dir, ifo, year, month, day)
                
                if verbose:
                    if os.path.exists(f'{dest_dir}/cwb.lcf'):
                        print(f"-> Done (cwb.lcf file)")
                    else:
                        raise IOError(f"The cwb.lcf file was NOT created")
                        
                if verbose:
                    if os.path.exists(f'{dest_dir}/EVENTS.txt'):
                        print(f"-> Done (EVENTS.txt file)")
                    else:
                        raise IOError(f"The EVENTS.txt file was NOT created")


                ''' Write Config '''
                
                if verbose:
                    print(f"\n|| Creating cwb_daily.ini file...")
                
                # Runs the Write_config function
                Write_config(dest_dir, config_file, ifo, UPV_flag_name, UPV_flag_file, hveto_flag_name, hveto_flag_file, paddings, UPV_channel_names, hveto_channel_names)

                if verbose:
                    if os.path.exists(f'{dest_dir}/cwb_daily.ini'):
                        print(f"-> Done (cwb_daily.ini file)")
                    else:
                        raise IOError(f"The cwb_daily.ini file was NOT created")


                ''' Make Executable '''

                if verbose:
                    print(f"\n|| Creating run.sh file...")

                # Runs the make_excecutable function
                make_excecutable(dest_dir, exe_file, config_file, year, month, day)
                
                if verbose:
                    if os.path.exists(f'{dest_dir}/run.sh'):
                        print(f"-> Done (run.sh file)")
                    else:
                        raise IOError(f"The run.sh file was NOT created")


                ''' Execute VET '''


                if verbose:
                    print(f"\n|| Executing VET...")

                # The command that executes VET
                cmd = 'cd {Dest_dir} && chmod +755 {Exe_file} && ./{Exe_file} '.format(Dest_dir=dest_dir, Exe_file=exe_file)
                os.system(cmd)

                if verbose:
                    print(f"-> Done (executing VET)")

                print(f'\nVET done with {year}, {month}, {day}\n')
                #print('type followings ...')
                #print('$ cd {Dest_dir} && chmod +755 {Exe_file} && ./{Exe_file} '.format(Dest_dir=dest_dir, Exe_file=exe_file))

        else:
            print(f"\nThere were no results on {year}, {month}, {day}\n")
            
        #no_data_counter += 1



    # if we are running over many days in a month.....
    else:
        for j in range(len(day)):

            # Converts the year, month, day to GPS date
            gps = to_gps(datetime.datetime(year, month, day[j]))
            date = gps.gpsSeconds
            date_end = date + 86400
            
            if verbose:
                print(f"Processing date: {year}-{month}-{day[j]}\n")
           
            if verbose:
                print(f"GPS time: {date}-{date_end}\n")


            # Path to the UPV veto segment round winners page
            directory = f"/home/detchar/public_html/NewUPV/O3/results/O3-H-{date}_{date_end}/"

            # Path for the html page
            file = directory + "index.html"


            if os.path.exists(file):

                #sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
                # custom modules
                #from auto_vet.utils.utils import query_hevto_flag, make_flag_xml, query_cwb_events, Write_config, make_excecutable # the list of glitch classes in each IFO

                work_dir = os.getcwd()

                if True:#__name__ == '__main__':

                    # make paddings
                    if pads is not None:
                        paddings = [[int(pad.split(',')[0][1:]), int(pad.split(',')[-1][:-1])] for pad in pads]
                    else:
                        paddings = None

                    # Converts the date to the correct format
                    date = "{}{:02d}{:02d}".format(year, month, day[j])

                    # Checks if the Flag Files already exist for UPV and hveto

                    # --------------------------------------------------------+
                    # Check if the {ifo}-UPV_ALLROUNDS.xml file already exits    
                    # --------------------------------------------------------+
                    if os.path.exists(f'{dest_dir}{ifo}-UPV_ALLROUNDS.xml'):
                        # If it already exits, remove it
                        os.remove(f'{dest_dir}{ifo}-UPV_ALLROUNDS.xml') 

                    # --------------------------------------------------------+
                    # Check if the {ifo}-HVETO_ALLROUNDS.xml file already exits    
                    # --------------------------------------------------------+
                    if os.path.exists(f'{dest_dir}{ifo}-HVETO_ALLROUNDS.xml'):
                        # If it already exits, remove it
                        os.remove(f'{dest_dir}{ifo}-HVETO_ALLROUNDS.xml') 


                    # make the destination directory if it does not exist yet
                    try:
                        os.mkdir(os.path.join(dest_dir))
                    except OSError:
                        pass


                    ''' UPV '''


                    # Runs the UPV Scraper function
                    channels, UPV_channel_names, UPV_flag_name, UPV_flag_file = UPV_scraper(year, month, day[j], dest_dir, ifo)    

                    if verbose:
                        print(f"|| Creating the UPV_sorted.txt file for all data...")

                    # Runs the UPV_sort_on_col function for the entire data
                    ascended_UPV = UPV_sort_on_col(data = np.loadtxt(f'{dest_dir}/UPV.txt', unpack = True), col = 1, ascending = True)
                    # Arranges the UPV vetoes in ascending order
                    start_time_ascend = ascended_UPV[0]
                    end_time_ascend = ascended_UPV[1]
                    # Creates the UPV_sorted.txt file containing the veto times in ascending order
                    np.savetxt(dest_dir + '/UPV_sorted_rd0.txt', np.transpose([start_time_ascend, end_time_ascend]), delimiter = ' ')

                    if verbose:
                        if os.path.exists(f'{dest_dir}/UPV_sorted_rd0.txt'):
                            print(f"-> Done (UPV_sorted_rd0.txt)")
                        else:
                            raise IOError("UPV_sorted_rd0.txt was NOT created")

                    if verbose:
                        print(f"|| Creating the UPV_sorted.txt file for each round...")

                    # Runs the UPV_sort_on_col function for every round
                    for i in range(1, channels + 1):
                        ascended_UPV = UPV_sort_on_col(data = np.loadtxt(f'{dest_dir}/UPV_rd{i}.txt', unpack = True), col = 1, ascending = True)
                        # Arranges the UPV vetoes in ascending order
                        start_time_ascend = ascended_UPV[0]
                        end_time_ascend = ascended_UPV[1]
                        # Creates the UPV_sorted.txt file containing the veto times in ascending order
                        np.savetxt(dest_dir + f'/UPV_sorted_rd{i}.txt', np.transpose([start_time_ascend, end_time_ascend]), delimiter = ' ')
                        
                        if verbose:
                            if os.path.exists(f'{dest_dir}/UPV_sorted_rd{i}.txt'):
                                print(f"-- Created UPV_sorted file for round {i} - {UPV_channel_names[i-1]}")
                            else:
                                raise IOError(f"UPV_sorted.txt for round {i} (UPV_sorted_rd{i}.txt) was NOT created")


                    ''' UPV Query Flag '''


                    # Runs the query_UPV_flag function
                    query_UPV_flag(year, month, day[j], dest_dir, ifo)  


                    ''' Hveto '''


                    hveto_channels, hveto_channel_names, hveto_flag_name, hveto_flag_file = hveto(year, month, day[j], dest_dir, ifo)

                    if verbose:
                        print(f"|| Creating the hveto_sorted.txt file for all data...")

                    # Runs the hveto_sort_on_col function for the entire data
                    ascended_hveto = hveto_sort_on_col(data = np.loadtxt(f'{dest_dir}/hveto.txt', unpack = True), col = 1, ascending = True)
                    # Arranges the hveto vetoes in ascending order
                    hveto_start_time_ascend = ascended_hveto[0]
                    hveto_end_time_ascend = ascended_hveto[1]
                    # Creates the hveto_sorted.txt file containing the veto times in ascending order
                    np.savetxt(dest_dir + '/hveto_sorted_rd0.txt', np.transpose([hveto_start_time_ascend, hveto_end_time_ascend]), delimiter = ' ')
                    
                    if verbose:
                        if os.path.exists(f'{dest_dir}/hveto_sorted_rd0.txt'):
                            print(f"-> Done (hveto_sorted_rd0.txt)")
                        else:
                            raise IOError("hveto_sorted_rd0.txt was NOT created")

                    if verbose:
                        print(f"|| Creating the hveto_sorted.txt file for each round...")


                    # Runs the hveto_sort_on_col function for every round
                    for k in range(1, hveto_channels + 1):
                        ascended_hveto = hveto_sort_on_col(data = np.loadtxt(f'{dest_dir}/hveto_rd{k}.txt', unpack = True), col = 1, ascending = True)
                        # Arranges the hveto vetoes in ascending order
                        hveto_start_time_ascend = ascended_hveto[0]
                        hveto_end_time_ascend = ascended_hveto[1]
                        # Creates the hveto_sorted.txt file containing the veto times in ascending order
                        np.savetxt(dest_dir + f'/hveto_sorted_rd{k}.txt', np.transpose([hveto_start_time_ascend, hveto_end_time_ascend]), delimiter = ' ')
                        
                        if verbose:
                            if os.path.exists(f'{dest_dir}/hveto_sorted_rd{k}.txt'):
                                print(f"-- Created hveto_sorted file for round {k} - {hveto_channel_names[k-1]}")
                            else:
                                raise IOError(f"hveto_sorted.txt for round {k} (hveto_sorted_rd{k}.txt) was NOT created")
      
                            
                    ''' Hveto Query Flag '''


                    # Runs the query_hveto_flag function
                    query_hveto_flag(year, month, day[j], dest_dir, ifo)



                    ''' Make Flag XML '''


                    if verbose:
                        print(f"\n|| Creating flag xml files for UPV...\n")

                    # Runs the make_flag_xml function for UPV - for all of the data and all rounds
                    for n in UPV_flag_name:
                        make_flag_xml(dest_dir, ifo, year, month, day[j], flag_name = UPV_flag_name[n], flag_file = UPV_flag_file[n], active_file = f'UPV_sorted_rd{n}.txt')
                        
                        if verbose:
                            if n == 0:
                                if os.path.exists(f'{dest_dir}/{UPV_flag_file[n]}'):
                                    print(f"\t-- xml file for all data created - {UPV_flag_file[n]}\n")
                                else:
                                    raise IOError(f"xml file for all data ({UPV_flag_file[n]}) was NOT created")
                            else:
                                if os.path.exists(f'{dest_dir}/{UPV_flag_file[n]}'):
                                    print(f"\t-- xml file for round {n} created - {UPV_flag_file[n]}\n")
                                else:
                                    raise IOError(f"xml file for round {n} ({UPV_flag_file[n]}) was NOT created")

                                    
                    if verbose:
                        print(f"-> Done (flag xml files for all rounds of UPV)")


                    if verbose:
                        print(f"\n|| Creating flag xml files for hveto...\n")

                    # Runs the make_flag_xml function for hveto - for all of the data and all rounds
                    for m in hveto_flag_name:
                        make_flag_xml(dest_dir, ifo, year, month, day[j], flag_name = hveto_flag_name[m], flag_file = hveto_flag_file[m], active_file = f'hveto_sorted_rd{m}.txt')
                        
                        if verbose:
                            if m == 0:
                                if os.path.exists(f'{dest_dir}/{hveto_flag_file[m]}'):
                                    print(f"\t-- xml file for all data created - {hveto_flag_file[m]}\n")
                                else:
                                    raise IOError(f"xml file for all data ({hveto_flag_file[m]}) was NOT created")
                            else:
                                if os.path.exists(f'{dest_dir}/{hveto_flag_file[m]}'):
                                    print(f"\t-- xml file for round {m} created - {hveto_flag_file[m]}\n")
                                else:
                                    raise IOError(f"xml file for round {m} ({hveto_flag_file[m]}) was NOT created")


                    if verbose:
                        print(f"-> Done (flag xml files for all rounds of hveto)")


                    ''' Query cWB Events '''

                    if verbose:
                        print(f"\n|| Creating cwb.lcf and EVENTS.txt files...")

                    # Runs the query_cwb_events function
                    query_cwb_events(dest_dir, ifo, year, month, day[j])
                    # print(dest_dir, ifo, year, month, day[j])
                    
                    if verbose:
                        if os.path.exists(f'{dest_dir}/cwb.lcf'):
                            print(f"-> Done (cwb.lcf file)")
                        else:
                            raise IOError(f"The cwb.lcf file was NOT created")
                            
                    if verbose:
                        if os.path.exists(f'{dest_dir}/EVENTS.txt'):
                            print(f"-> Done (EVENTS.txt file)")
                        else:
                            raise IOError(f"The EVENTS.txt file was NOT created")


                    ''' Write Config '''

                    if verbose:
                        print(f"\n|| Creating cwb_daily.ini file...")

                    # Runs the Write_config function
                    Write_config(dest_dir, config_file, ifo, UPV_flag_name, UPV_flag_file, hveto_flag_name, hveto_flag_file, paddings, UPV_channel_names, hveto_channel_names)
                    
                    if verbose:
                        if os.path.exists(f'{dest_dir}/cwb_daily.ini'):
                            print(f"-> Done (cwb_daily.ini file)")
                        else:
                            raise IOError(f"The cwb_daily.ini file was NOT created")


                    ''' Make Executable '''

                    if verbose:
                        print(f"\n|| Creating run.sh file...")

                    # Runs the make_excecutable function
                    make_excecutable(dest_dir, exe_file, config_file, year, month, day[j])

                    if verbose:
                        if os.path.exists(f'{dest_dir}/run.sh'):
                            print(f"-> Done (run.sh file)")
                        else:
                            raise IOError(f"The run.sh file was NOT created")
                    

                    ''' Execute VET '''


                    if verbose:
                        print(f"\n|| Executing VET...")

                    # The command that executes VET
                    cmd = 'cd {Dest_dir} && chmod +755 {Exe_file} && ./{Exe_file} '.format(Dest_dir=dest_dir, Exe_file=exe_file)
                    os.system(cmd)

                    if verbose:
                        print(f"-> Done (executing VET)")

                    print(f'\n**VET done with {year}, {month}, {day[j]}**\n')
                    #print('type followings ...')
                    #print('$ cd {Dest_dir} && chmod +755 {Exe_file} && ./{Exe_file} '.format(Dest_dir=dest_dir, Exe_file=exe_file))

            else:
                print(f"\nThere were no results on {year}, {month}, {day[j]}\n")
                timer = True
                #no_data_counter += 1

            # After a day is finished, removes all files created to make sure directories don't get overfilled
            os.system(f'find {dest_dir} -maxdepth 1 -type f -delete')



    timer_end = time.time()
    total_time = (timer_end - timer_start) / 60
    print(f"Total run time: {total_time:.2f} minutes")
        
    #print(f"There were {no_data_counter} days with no data")


''' End of Wrapper() '''

#=================================================================================================================
#=================================================================================================================