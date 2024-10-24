import skrf as rf
from skrf.media.cpw import CPW

def assemble_network(options):
    network_type = options["type"]
    if network_type == "simple":
        return simple_network(options)
    elif network_type == "feedline_resonator":
        return tl_resonator_network(options)
    elif network_type == "two_coupled":
        return two_coupled_qubits_network(options)


def simple_network(options):
    f_min, f_max = options["frequency"]
    freq = rf.Frequency(start=f_min, stop=f_max, unit='GHz', npoints=options["n_points"])
    capacitor = rf.DefinedGammaZ0(freq, z0=50).capacitor
    inductor = rf.DefinedGammaZ0(freq, z0=50).inductor
    a, b, ep_r = options["a"], options["b"], options["ep_r"]
    cpw = CPW(frequency = freq, w = a, s = b, ep_r = ep_r)
    port1 = rf.Circuit.Port(freq, 'port1', z0=50) # Launcher 1
    port2 = rf.Circuit.Port(freq, 'port2', z0=50) # Launcher 2

    tl = []
    cplr = []
    qplr = []

    readout_res_c = []
    readout_res_l = []
    qubit_c = []
    qubit_l = []
    gnd = []

    cnx = []
    
    for i, fl in enumerate(options["feedline_traits"]):
        tl.append(cpw.line(d = fl/1e6, unit='m', name='trans_line_'+str(i)))
    for i, (ck, cqr, rr_c, rr_l, qb_c, qb_l) in enumerate(zip(
                                                                options["ck"],
                                                                options["cqr"],                                                  
                                                                options["readout_resonator_c"],
                                                                options["readout_resonator_l"],
                                                                options["qubit_c"],
                                                                options["qubit_l"])):
        cplr.append(capacitor(ck, name='ck_'+str(i)))
        qplr.append(capacitor(cqr, name='cqr_'+str(i)))
        readout_res_c.append(capacitor(rr_c, name='rr_c_'+str(i)))
        readout_res_l.append(inductor(rr_l, name='rr_l_'+str(i)))
        qubit_c.append(capacitor(qb_c, name='qb_c_'+str(i)))
        qubit_l.append(inductor(qb_l, name='qb_l_'+str(i)))
    n_qubits = len(options["qubit_c"])    
    for i in range(len(options["qubit_c"]) + 1):
        gnd.append(rf.Circuit.Ground(freq, name='gnd_'+str(i)))

    n_qubits = len(options["qubit_c"])
    for i in range(n_qubits+1):
        if i == 0:
            cnx.append([(port1, 0), (tl[i], 0)])
        elif i == n_qubits:
            cnx.append([(tl[i], 1), (port2, 0)])
            break
        # else:
        cnx.append([(tl[i], 1), (tl[i+1], 0), (cplr[i], 0),])
        cnx.append([(cplr[i], 1), (readout_res_c[i], 0), (readout_res_l[i], 0), (qplr[i], 0)])
        cnx.append([(qplr[i], 1), (qubit_c[i], 0), (qubit_l[i], 0)])
        cnx.append([(readout_res_c[i], 1), (readout_res_l[i], 1), (qubit_c[i], 1), (qubit_l[i], 1), (gnd[i], 0)])

    return cnx

def tl_resonator_network(options):
    # In progress
    f_min, f_max = options["frequency"]
    freq = rf.Frequency(start=f_min, stop=f_max, unit='GHz', npoints=options["n_points"])
    capacitor = rf.DefinedGammaZ0(freq, z0=50).capacitor
    inductor = rf.DefinedGammaZ0(freq, z0=50).inductor
    a, b, ep_r = options["a"], options["b"], options["ep_r"]
    cpw = CPW(frequency = freq, w = a, s = b, ep_r = ep_r)
    port1 = rf.Circuit.Port(freq, 'port1', z0=50) # Launcher 1
    port2 = rf.Circuit.Port(freq, 'port2', z0=50) # Launcher 2

    tl = []
    cplr = []
    qplr = []

    readout_res_c = []
    readout_res_l = []
    qubit_c = []
    qubit_l = []
    gnd = []    
    for i, fl in enumerate(options["feedline_traits"]):
        tl.append(cpw.line(d = fl/1e6, unit='m', name='trans_line_'+str(i)))
    end_tl = cpw.line(d = 500/1e6, unit='m', name='end_tl')
    for i, (ck, cqr, rr_c, rr_l, qb_c, qb_l) in enumerate(zip(
                                                                options["ck"],
                                                                options["cqr"],                                                  
                                                                options["readout_resonator_c"],
                                                                options["readout_resonator_l"],
                                                                options["qubit_c"],
                                                                options["qubit_l"])):
        cplr.append(capacitor(ck, name='ck_'+str(i)))
        qplr.append(capacitor(cqr, name='cqr_'+str(i)))
        readout_res_c.append(capacitor(rr_c, name='rr_c_'+str(i)))
        readout_res_l.append(inductor(rr_l, name='rr_l_'+str(i)))
        qubit_c.append(capacitor(qb_c, name='qb_c_'+str(i)))
        qubit_l.append(inductor(qb_l, name='qb_l_'+str(i)))
        
    feedline_capacitor = capacitor(options["cfl"], name = "cfl")
    n_qubits = len(options["qubit_c"])    
    for i in range(len(options["qubit_c"]) + 1):
        gnd.append(rf.Circuit.Ground(freq, name='gnd_'+str(i)))

    n_qubits = len(options["qubit_c"])
    cnx = []
    cnx.append([(port1, 0), (tl[0], 0)])
    cnx.append([(tl[0], 1), (feedline_capacitor, 0)])
    for i in range(n_qubits+1):
        if i == 0:
            cnx.append([(feedline_capacitor, 1), (tl[i+1], 0)])
        elif i == n_qubits:
            cnx.append([(tl[i+1], 1), (tl[i+2], 0), (tl[i+3], 0),])
            cnx.append([(tl[i+2], 1), (port2, 0)])
            cnx.append([(tl[i+3], 1), (gnd[i], 0)])
            break
        # else:
        cnx.append([(tl[i+1], 1), (tl[i+2], 0), (cplr[i], 0),])
        cnx.append([(cplr[i], 1), (readout_res_c[i], 0), (readout_res_l[i], 0), (qplr[i], 0)])
        cnx.append([(qplr[i], 1), (qubit_c[i], 0), (qubit_l[i], 0)])
        cnx.append([(readout_res_c[i], 1), (readout_res_l[i], 1), (qubit_c[i], 1), (qubit_l[i], 1), (gnd[i], 0)])      
    return cnx

def two_coupled_qubits_network(options):
    f_min, f_max = options["frequency"]
    freq = rf.Frequency(start=f_min, stop=f_max, unit='GHz', npoints=options["n_points"])
    capacitor = rf.DefinedGammaZ0(freq, z0=50).capacitor
    inductor = rf.DefinedGammaZ0(freq, z0=50).inductor
    a, b, ep_r = options["a"], options["b"], options["ep_r"]
    cpw = CPW(frequency = freq, w = a, s = b, ep_r = ep_r)
    port1 = rf.Circuit.Port(freq, 'port1', z0=50) # Launcher 1
    port2 = rf.Circuit.Port(freq, 'port2', z0=50) # Launcher 2

    tl = []
    cplr = []
    qplr = []

    readout_res_c = []
    readout_res_l = []
    qubit_c = []
    qubit_l = []
    gnd = []
    for i in range(len(options["qubit_c"]) + 1):
        gnd.append(rf.Circuit.Ground(freq, name='gnd_'+str(i)))

    cnx = []
    
    for i, fl in enumerate(options["feedline_traits"]):
        tl.append(cpw.line(d = fl/1e6, unit='m', name='trans_line_'+str(i)))
    for i, (ck, cqr, rr_c, rr_l, qb_c, qb_l) in enumerate(zip(
                                                                options["ck"],
                                                                options["cqr"],                                                  
                                                                options["readout_resonator_c"],
                                                                options["readout_resonator_l"],
                                                                options["qubit_c"],
                                                                options["qubit_l"])):
        cplr.append(capacitor(ck, name='ck_'+str(i)))
        qplr.append(capacitor(cqr, name='cqr_'+str(i)))
        readout_res_c.append(capacitor(rr_c, name='rr_c_'+str(i)))
        readout_res_l.append(inductor(rr_l, name='rr_l_'+str(i)))
        qubit_c.append(capacitor(qb_c, name='qb_c_'+str(i)))
        qubit_l.append(inductor(qb_l, name='qb_l_'+str(i)))

    ccplr_0 = capacitor(options["cc"][0], name = 'cc_0')
    ccplr_1 = capacitor(options["cc"][1], name = 'cc_1')
    coupler_res_c = capacitor(options["coupler_res_c"], name = 'cp_res_c')
    coupler_res_l = inductor(options["coupler_res_l"], name = 'cp_res_l')
    # gnd.append(rf.Circuit.Ground(freq, name='gnd_coupler'))

    cnx = [
            [(port1, 0), (tl[0], 0)],
            [(tl[0], 1), (tl[1], 0), (cplr[0], 0)],
            [(tl[1], 1), (tl[2], 0), (cplr[1], 0)],
            [(tl[2], 1), (port2, 0)],



            [(cplr[0], 1), (readout_res_c[0], 0), (readout_res_l[0], 0), (qplr[0], 0)], 
            [(qplr[0], 1), (qubit_c[0], 0), (qubit_l[0], 0), (ccplr_0, 0)],
            [(readout_res_c[0], 1), (readout_res_l[0], 1), (qubit_c[0], 1), (qubit_l[0], 1), (gnd[0], 0)],

            [(ccplr_0, 1), (coupler_res_c, 0), (coupler_res_l, 0), (ccplr_1, 1)],
            [(coupler_res_c, 1), (coupler_res_l, 1), (gnd[2], 0)],

            [(cplr[1], 1), (readout_res_c[1], 0), (readout_res_l[1], 0), (qplr[1], 0)], 
            [(qplr[1], 1), (qubit_c[1], 0), (qubit_l[1], 0), (ccplr_1, 0)],
            [(readout_res_c[1], 1), (readout_res_l[1], 1), (qubit_c[1], 1), (qubit_l[1], 1), (gnd[1], 0)],
            ]


    return cnx