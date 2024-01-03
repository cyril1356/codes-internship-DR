
def give_electron_config_stars(closed_config, transition, end_stars):
    """
    gives all stars configurastion
    transition: explain the transition on witch the bound electron will go ex 5p
    list_states ground: are all the values configuration for n_ground: exemple for the transition 5p we have to put
    end star is an int saying the shell in witch the d electron will be catch
    """

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g"
    config = config.split()
    closed_config = closed_config.split()
    indice_bound_transition = config.index(transition)
    value_n_bound = config[indice_bound_transition][0]
    list_bound = []

    # bound electrons
    for element in config:
        if element not in closed_config:
            if element[0] == value_n_bound:
                if element == config[indice_bound_transition]:
                    list_bound.append(transition + "2 ")
                else:
                    if config.index(element) > indice_bound_transition:
                        list_bound.append(
                            config[indice_bound_transition] + " " + element
                        )
                    else:
                        list_bound.append(
                            element + " " + config[indice_bound_transition]
                        )

    for i in range(int(value_n_bound) + 1, end_stars + 1):
        list_bound.append(str(i) + "* ")
    return list_bound


print(give_electron_config_stars("1s 2s 2p 3s 3p 4s 3d 4p", "5p", 8))