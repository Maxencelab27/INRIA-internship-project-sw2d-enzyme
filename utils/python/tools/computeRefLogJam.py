import csv
import classTransportType

def generateur_timer(begin, end, range):
    """
    Generate an iterator on time, format [d, h, min, sec, msec]
    """
    timer = begin
    while timer <= end:
        yield(timer)

        test = [elm1 + elm2 for (elm1, elm2) in zip(timer, range)]
        if test[4] >= 1000:
            test[4] = test[4] % 1000
            test[3] += 1
        if test[3] >= 60:
            test[3] = test[3] % 60
            test[2] += 1
        if test[2] >= 60:
            test[2] = test[2] % 60
            test[1] += 1
        if test[1] >= 24:
            test[1] = test[1] % 24
            test[0] += 1
        timer = test

def convert_time_iso8601(time):
    """
    Convert a time [d, h, m, sec, msec] onto the string format iso8601
    """
    nbr_days = str(time[0])
    nbr_days = '0'*(4-len(nbr_days)) + nbr_days
    nbr_heure = str(time[1])
    nbr_heure = '0'*(2-len(nbr_heure)) + nbr_heure
    nbr_min = str(time[2])
    nbr_min = '0'*(2-len(nbr_min)) + nbr_min
    nbr_sec = str(time[3])
    nbr_sec = '0'*(2-len(nbr_sec)) + nbr_sec
    nbr_msec = str(time[4])
    nbr_msec = '0'*(3-len(nbr_msec)) + nbr_msec

    return nbr_days + 'd' + nbr_heure+'h'+nbr_min+'m'+nbr_sec+'s'+nbr_msec

def write_header(output_filename):
    """
    Write LogJam header, using the same format than sw2d
    """
    with open(output_filename, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t')
        spamwriter.writerow(["### LogJam ###"])
        spamwriter.writerow(["Id", "x", "y", "v_x", "v_y", "theta", "omega", "rho", "length"])

def write_cylinder(output_filename, my_cylinder : classTransportType.Cylinder):
    """
    Write a Cylinder in output_filename, using only its last position/velocities    
    """
    with open(output_filename, 'a', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t')
        spamwriter.writerow(["%i"%(my_cylinder.id),
                            "%.5f"%(my_cylinder.position[-1][0]),
                            "%.5f"%(my_cylinder.position[-1][1]),
                            "%.5f"%(my_cylinder.velocity[-1][0]),
                            "%.5f"%(my_cylinder.velocity[-1][1]),
                            "%.5f"%(my_cylinder.theta[-1]),
                            "%.5f"%(my_cylinder.omega[-1]),
                            "%.5f"%(my_cylinder.rho),
                            "%.5f"%(my_cylinder.length)])
