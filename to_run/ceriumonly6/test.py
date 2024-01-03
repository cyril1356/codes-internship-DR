

def give_electron_config_FAC(closed_config: str, d_state: str):
    """
    function made for fac calculations:
    wath out closed config no comma only spaces, closed config put the 'i' state inside too.

    exemple:
    give_electron_config_FAC('1s 2s 2p 3s 3p 4s 3d 4p 5s','5p 5g') ==> ['4d 5g', '5p 8s', '5p 7p', '5p 6d', '5p 5f', '5p 7s', '5p 6p', '5p 5d', '5p 4f', '5p 6s', '5p 5p']
    """

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"
    config = config.replace(closed_config, "")
    config = config.split()
    d_state = d_state.split()
    ##bound electron
    bound = d_state[0]
    indice_stop1 = config.index(bound)
    indice_ini = config.index(d_state[1])
    res = []
    for i in range(indice_stop1):
        a = config[i] + " " + d_state[1]
        res += [a]
    # free electron

    for i in range(
        indice_ini - 1, indice_stop1, -1
    ):  ##put -1 if you want have config 5p 5p for d_state[0]=5p
        res += [d_state[0] + " " + config[i]]
    return res




print(give_electron_config_FAC("1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p","6p 5g"))