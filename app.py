from flask import Flask, render_template, request
from Bio.Blast import NCIWWW, NCBIXML
from Bio import SeqIO
import smtplib
import ssl
from email.message import EmailMessage
import OS
import csv
import tempfile
app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def blast():
    if request.method == "POST":
        #get the user's email
        email = request.form["email"]

        #get pasted FASTA sqeuences
        fasta_text = request.form.get("fasta","").strip()

        #get file
        fasta_file = request.files.get("fasta_file")

        #combine all input sequences into a single list
        sequences = []
        if fasta_text:
            sequences.extend(list(SeqIO.parse(fasta_text.splitlines(), "fasta")))
        if fasta_file and fasta_file.filename:
            sequences.extend(list(SeqIO.parse(fasta_file.stream, "fasta")))
        if not sequences:
            return "No valid sequences given!"
        
        #creating temporary csv file to write the blast results
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w", newline="") as output_csv:
            writer = csv.writer(output_csv)
            writer.writerow([
                "Query ID", " Hit ID", "Hit Def", "Identity %", "Coverage %", "Gaps", "Alginment Length", "E-value", "Sequence"
            ])

            #loop over each protein sequence
            for record in sequences:
                #send to BLAST remotely against nr
                result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
                blast_record = NCBIXML.read(result_handle)
                query_len = blast_record.query_length

                #extract data from BLAST RESULTS
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        identity_pct = round((hsp.identities / hsp.align_length )*100, 2)
                        coverage_pct = round((hsp.align_length / query_len )*100, 2)
                        writer.writerow([
                            blast_record.query_id,
                            alignment.hit_id,
                            alignment.hit_def,
                            identity_pct,
                            coverage_pct,
                            hsp.gaps,
                            hsp.align_length,
                            hsp.expect,
                            hsp.sbjct  #aligned sequences
                        ])
        #email the results to the user
        send_email_with_attachment(
            sender="default@gmail.com",
            recipient = email,
            filepath=output_csv.name
        )

        #delete temp file
        os.remove(output_csv.name)

        return render_template("result.html", email=email)
    #render the upload form
    return render_template("index.html")

def send_email_with_attachment(sender, recipient, filepath):
    # Replace with a real mail and password
    password = "password"

    #prepare email
    msg = EmailMessage()
    msg["Subject"] = "Your BLAST Results"
    msg["From"] = sender
    msg["To"]= recipient
    msg.set_content("Attached are your BLAST RESULTS in CSV format")

    #attach csv file
    with open(filepath, "rb") as f:
        data = f.read()
        msg.add_attachment(data, maintype= "application", subtype="csv", filename=os.path.basename(filepath))

        #securely connect to GMAIL SMTP and sen the mail
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gamil.com", 465, context=context) as smtp:
            smtp.login(sender, password)
            smtp.send_message(msg)
            

