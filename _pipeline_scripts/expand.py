
def Expand(Rr, Br, Al, At):
    for r in Rr:
        for b in Br:
            for l in Al:
                for t in At:
                    print "%s\t%s\t%s\t%s\n"  % (r,b,l,t)


Expand(["Rr"], ["Br"], ["Al"], ["At"])
Expand(["Rr"], ["Br"], ["AlA", "AlB"], ["AtA", "AtB"])

