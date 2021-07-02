import json


data = json.load(open("./reference/sars_cov_2.json"))

pfam_domains = []
for g in data["genes"]:
    for t in g["transcripts"]:
        for t2 in t["translations"]:
            for p in t2["protein_features"]:
                if p["dbname"] == "Pfam":
                    pfam_domains.append({
                        "start": int(g["start"]) + int(p["start"]) - 1,     # from 1-based to 0-based
                        "end": int(g["start"]) + int(p["end"]) - 1,
                        "interpro_name":p["interpro_name"],
                        "interpro_description": p["interpro_description"]
                    })

chromosome = "MN908947.3"

with open("reference/pfam_names.bed", "w") as f:
    for d in pfam_domains:
        f.write("\t".join([chromosome, str(d["start"]), str(d["end"]), d["interpro_name"]]) + "\n")

with open("reference/pfam_descriptions.bed", "w") as f:
    for d in pfam_domains:
        f.write("\t".join([chromosome, str(d["start"]), str(d["end"]), d["interpro_description"]]) + "\n")
