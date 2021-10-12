import xml.etree.cElementTree as ET

import numpy as np


def parse_varray(varray):
    if varray.get("type") == "int":
        m = np.array(
            [[int(number) for number in v.text.split()] for v in varray.findall("v")],
            dtype=np.int64,
        )
    else:
        m = np.array(
            [[float(number) for number in v.text.split()] for v in varray.findall("v")],
            dtype=np.float64,
        )
    return m


def _parse_calculation(elem):
    try:
        istep = {
            i.attrib["name"]: float(i.text) for i in elem.find("energy").findall("i")
        }
    except AttributeError:  # not all calculations have an energy
        istep = {}
        pass

    struct = elem.find("structure")
    pos = parse_varray(struct.find("varray"))  # direct
    lat = parse_varray(struct.find("crystal").find("varray"))
    # pos = np.dot(pos, lat) # direct to cart

    for va in elem.findall("varray"):
        istep[va.attrib["name"]] = parse_varray(va)

    istep["positions"] = pos
    istep["lat"] = lat
    return istep


def parse_incar(fn):
    params = {}
    ee = ET.iterparse(fn)
    for _, elem in ee:
        tag = elem.tag
        if tag == "incar":
            for c in elem:
                name = c.attrib.get("name")
                val = c.text.strip() if c.text else ""
                params[name] = val
            break
    return params


def parse_force_pos(fn, skip, nsteps):
    ee = ET.iterparse(fn)
    count = 0
    count_read = 0

    try:
        for _, elem in ee:
            tag = elem.tag
            if tag == "calculation":
                count += 1
                if count >= skip and count_read < nsteps:
                    cal = _parse_calculation(elem)
                    count_read += 1
                    yield (count, count_read, cal["forces"], cal["positions"])
                elif count_read >= nsteps:
                    raise StopIteration
    except ET.ParseError:
        print("XML is malformed. Parsing has stopped but partial data is available.")
        raise StopIteration


def count_sc_steps(fn):
    ee = ET.iterparse(fn)
    res = []
    try:
        for _, elem in ee:
            tag = elem.tag
            if tag == "calculation":
                n_sc_steps = len(elem.findall("scstep"))
                res.append(n_sc_steps)
        return res
    except ET.ParseError:
        # Calculation not success
        return None


def write(dat, fn):
    with open(fn, "ab") as f1:
        np.savetxt(f1, dat, fmt="%12.8f")

