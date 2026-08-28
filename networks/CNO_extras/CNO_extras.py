import pynucastro as pyna


def create_network():

    net = pyna.network_helper(["h1", "he4",
                               "c12", "c13",
                               "n13", "n14", "n15",
                               "o14", "o15", "o16", "o17", "o18",
                               "f17", "f18", "f19",
                               "ne18", "ne19", "ne20",
                               "mg22", "mg24"],
                              tabular_ordering=["ffn", "langanke", "oda"],
                              inert_nuclei=["fe56"], network_type="amrex")

    # remove the C-burning rates
    rr = net.get_rate_by_name(["c12(c12,a)ne20",
                               "o16(c12,a)mg24",
                               "ne20(a,c12)c12",
                               "mg24(a,c12)o16"])

    net.remove_rates(rr)

    # add additional breakout rates
    rl = pyna.ReacLibLibrary()

    rne18ag = rl.get_rate_by_name("ne18(a,g)mg22")
    net.remove_rates(rne18ag)

    rne18_mg22, _ = pyna.rates.make_ap_pg_rates(rl, "ne18", "mg22",
                                                use_detailed_balance=True)

    rne19pg = rl.get_rate_by_name("ne19(p,g)na20")
    rne19_mg22 = pyna.ModifiedRate(rne19pg,
                                   new_products=[pyna.Nucleus("mg22")],
                                   stoichiometry={pyna.Nucleus("p"): 3},
                                   description="Ne19(p,γ)Na20(p,γ)Mg21(β+)Na21(p,γ)Mg22")

    rne20pg = rl.get_rate_by_name("ne20(p,g)na21")
    rne20_mg22 = pyna.ModifiedRate(rne20pg,
                                   new_products=[pyna.Nucleus("mg22")],
                                   stoichiometry={pyna.Nucleus("p"): 2},
                                   description="Ne20(p,γ)Na21(p,γ)Mg22")

    net.add_rates([rne18_mg22, rne19_mg22, rne20_mg22])

    return net


def doit():

    net = create_network()

    net.write_network()

    comp = pyna.Composition(net.get_nuclei(), init="solar")

    rho = 1.e4
    T = 1.e8

    net.plot(rho, T, comp, outfile="cno_extras.png",
             node_size=500, node_font_size="9",
             Z_range=[1, 13], N_range=[1, 13])

    net.plot(outfile="cno_extras_hide_alpha.png",
             node_size=500, node_font_size="9",
             Z_range=[1, 13], N_range=[1, 13],
             rotated=True,
             hide_xalpha=True)

    net.summary()

if __name__ == "__main__":
    doit()
