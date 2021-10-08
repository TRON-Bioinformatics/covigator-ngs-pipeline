import json


data = json.load(open("./sars_cov_2.json"))

pfam_domains = []
for g in data["genes"]:
    for t in g["transcripts"]:
        for t2 in t["translations"]:
            for p in t2["protein_features"]:
                if p["dbname"] == "Pfam":
                    pfam_domains.append({
                        # domains coordinates are in the protein space
                        # from 1-based to 0-based
                        "start": int(g["start"]) + (int(p["start"]) * 3) - 1,
                        "end": int(g["start"]) + (int(p["end"]) * 3) - 1,
                        "interpro_name":p["interpro_name"],
                        "interpro_description": p["interpro_description"]
                    })

chromosome = "MN908947.3"

pfam_domains = sorted(pfam_domains, key=lambda x: x["start"])

with open("pfam_names.bed", "w") as f:
    for d in pfam_domains:
        f.write("\t".join([chromosome, str(d["start"]), str(d["end"]), d["interpro_name"]]) + "\n")

with open("pfam_descriptions.bed", "w") as f:
    for d in pfam_domains:
        f.write("\t".join([chromosome, str(d["start"]), str(d["end"]), d["interpro_description"]]) + "\n")
