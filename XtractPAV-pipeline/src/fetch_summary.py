#-------------------------------------------------------------------------------
# Script for curate the report of PAV analysis
# Authors:
#    Rana Sheraz Ahmad           : ranasheraz.202101902@gcuf.edu.pk
#    Mr. Muhammad Sadaqat        : muhammad.sadaqat@univ-rennes.fr
#    Dr. Muhammad Tahir ul Qamar : tahirulqamar@gcuf.edu.pk
# Created Time: 
#        Mon Jun 02 23:49:05 2025
# Version: 
#        1.0 [ initial release ]
#-------------------------------------------------------------------------------
import sys
import plotly.graph_objects as go
import plotly.express as px

def length(bed, gff3, summary):
    table_rows = []  
    pav_lengths = [] 
    chrom_counts = {}  
    genic_count = {}

    with open(gff3, 'r') as gene:
        for line in gene:
            if line.startswith("#"): 
                continue
            
            columns = line.strip().split('\t')
            if len(columns) < 3:  
                continue
            
            accession = columns[0] 
            feature_type = columns[2] 
            
            if feature_type.lower() == "gene":  
                genic_count[accession] = genic_count.get(accession, 0) + 1


                

    with open(bed, 'r') as coords:
        pav_counter = 1
        
        for i in coords:
            if i.startswith("#") or i.strip() == "":
                continue

            parts = i.strip().split('\t')
            chrom = parts[0]  
            start = int(parts[1])
            end = int(parts[2])
            pav_length = end - start + 1  
            pav_lengths.append(pav_length) 

           
            if chrom in chrom_counts:
                chrom_counts[chrom] += 1
            else:
                chrom_counts[chrom] = 1

     
            table_rows.append(f"<tr><td>PAV{pav_counter}</td><td>{chrom}</td><td>{start}</td><td>{end}</td><td>{pav_length}</td></tr>")
            
            pav_counter += 1  

    if pav_lengths:
        
        longest_pav = max(pav_lengths) 
        shortest_pav = min(pav_lengths)
        average_pav = sum(pav_lengths) / len(pav_lengths)
        chart_max_y = longest_pav * 2  
        total_pavs = len(pav_lengths)
    else:
        longest_pav = shortest_pav = average_pav = chart_max_y = 0


    fig1 = go.Figure()
    fig1.add_trace(go.Bar(
        x=[f"PAV {i+1}" for i in range(len(pav_lengths))],  
        y=pav_lengths,
        marker_color='#003153',
        name="PAV Lengths"
    ))

    fig1.update_layout(
        title="PAV Length Distribution along with the their number",
        xaxis_title="PAV Number",
        yaxis_title="PAV Length (bp)",
       yaxis=dict(range=[0, chart_max_y], showgrid=False),
        bargap=0.05,  
        bargroupgap=0.05,
        template="plotly_white",
        xaxis=dict(showgrid=False),


    )

    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        y=list(chrom_counts.keys()),  
        x=list(chrom_counts.values()),  
        orientation='h',
        marker_color='#003153',
        name="PAV Counts"
    ))
    fig2.add_trace(go.Bar(
        y=list(genic_count.keys()),  
        x=[genic_count.get(chrom, 0) for chrom in genic_count.keys()],  
        orientation='h',
        marker_color='#960018',
        name="Genic PAVs"
    ))

    fig2.update_layout(
    title="Number of PAVs Per Chromosome",
    xaxis_title="Number of PAVs",
    yaxis_title="Chromosome",
    barmode='group', 
    bargap=0.1,  
    bargroupgap=0.02,
    template="plotly_white",
    plot_bgcolor="white",  
    paper_bgcolor="white",  
    xaxis=dict(showgrid=False),  
    yaxis=dict(showgrid=False)   
)

    fig3 = go.Figure()

    fig3.add_trace(go.Pie(
        labels=list(chrom_counts.keys()),  
        values=list(chrom_counts.values()),
        textinfo='percent+label',  
        marker=dict(colors=px.colors.sequential.Inferno),  
        hole=0.4,  
        name="PAV Distribution"
                    ))

    fig3.update_layout(
        title="Percentage of PAVs distribution per Chromosome",
        template="plotly_white"
                    )


    fig1_json = fig1.to_json()
    fig2_json = fig2.to_json()
    fig3_json = fig3.to_json()

    chrom_table_rows = "".join(f"<tr><td>{chrom}</td><td>{count}</td></tr>" for chrom, count in chrom_counts.items())
    genic_rows = "".join(f"<tr><td>{accession}</td><td>{count}</td></tr>" for accession, count in genic_count.items())

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PAVs Analysis Summary</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{
                background-color: #2C3539; 
                color: #FFFFFF; 
                font-family: Arial, sans-serif;
                padding: 20px;
            }}
            table {{ 
                border-collapse: collapse; 
                width: 100%;
                background-color: #36454F; 
            }}
            th, td {{ 
                border: 1px solid #FFFFFF; 
                padding: 8px; 
                text-align: center; 
            }}
            th {{ 
                background-color: #003153; 
                color: #FFFFFF; 
                # #556B2F
            }}
            .chart-container {{ 
                width: 100%; 
                height: 500px; 
                background-color: #4F545A; 
                border-radius: 10px;
                padding: 10px;
            }}
        </style>
    </head>
    <body>
        <h2>PAVs Analysis Summary</h2>
        <p><strong>Total PAVs founded :</strong> {total_pavs}</p>
        <p><strong>The length of Longest PAV is :</strong> {longest_pav} bp</p>
        <p><strong>The length of Shortest PAV is :</strong> {shortest_pav} bp</p>
        <p><strong>The Average PAV Length:</strong> {average_pav:.2f} bp</p>
        
        <h3>PAV Length Visualization</h3>
        <div class="chart-container" id="pav_chart"></div>
        <script>
            var fig1 = {fig1_json};
            Plotly.newPlot('pav_chart', fig1.data, fig1.layout);
        </script>

        <h3>PAV Counts Per Chromosome</h3>
        <div class="chart-container" id="pav_counts_chart"></div>
        <script>
            var fig2 = {fig2_json};
            Plotly.newPlot('pav_counts_chart', fig2.data, fig2.layout);
        </script>
	
	<h3>Pie-Chart of PAV Distribution by Chromosome</h3>
        <div class="chart-container" id="pav_pie_chart"></div>
        <script>
            var fig3 = {fig3_json};
            Plotly.newPlot('pav_pie_chart', fig3.data, fig3.layout);
        </script>


        <h3>Number of PAV on their chromosome</h3>
        <table>
            <tr>
                <th>Chromosome</th>
                <th>Total PAVs</th>
            </tr>
            {chrom_table_rows}
        </table>

        <h3>Genic PAVs </h3>
         <table>
            <tr>
                <th>Chromosome/Scafold</th>
                <th>Total Genic PAVs</th>
            </tr>
            {genic_rows}
        </table>

        <h3>PAV Details</h3>
        <table>
            <tr>
                <th>Number of PAV</th>
                <th>Chromosome</th>
                <th>Start</th>
                <th>End</th>
                <th>PAV Length</th>
            </tr>
            {''.join(table_rows)}
        </table>
    </body>
    </html>
"""



    with open(summary, "w", encoding="utf-8") as html_file:
        html_file.write(html_content)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <bed_file> <summary_html>")
        sys.exit(1)

    bed_file = sys.argv[1]
    gff3_file = sys.argv[2]
    summary = sys.argv[3]
    length(bed_file, gff3_file, summary)
