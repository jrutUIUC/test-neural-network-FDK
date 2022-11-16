<p><span style="color:#000000"><span style="font-family:Verdana,Geneva,sans-serif"><span style="font-size:16px"><strong>test neural network FDK</strong></span></span></span></p>

<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">This repository contains scipts and data used to evaluate FDK phenotyping based on a trained neural network. We evaluated how data generated by the neural network would impact genomic selection accuracy and we estimated trait heritabilities and genetic correlations</span></span></span></p>

<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">Data and R scrips in this repository are used for the analyses described in the manusrctipt entitled &#39;A Neural Network for Phenotyping Fusarium Damaged Kernels (FDK) and its Impact on Genomic Selection Accuracy&#39; </span></span></span></p>

<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">The genotypic data, in the file &#39;FHBgeno_2020-21.csv&#39; are arranged with lines in rows and markers in columns. Marker names are in the headers. The phenotypic data, in the file &#39;FHBpheno_2020-21.xls&#39; are arranged with column headers in the first line, and each row corresponds to a single experimental unit. </span></span></span></p>

<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">The headers are defined as follows:</span></span></span></p>

<table cellpadding="0" cellspacing="0" class="Table" style="border-collapse:collapse; width:751px">
	<tbody>
		<tr>
			<td style="border-color:black; border-style:solid; border-width:1px; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif"><strong>header</strong></span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:1px solid black; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif"><strong>description</strong></span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">GS_trainingset</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">A binary variable indicating whether the row of data was included in the genomic selection (GS) model training set.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">GS_validset</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;A binary variable indicating whether the row of data was included in the GS validation set.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">observationUnitDbId</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The database ID for the experimental unit in the breeding program database.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">studyYear</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The year when data were collected.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">germplasmName</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The name of the line being evaluated.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">synonym</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;A synonym for the name of the line being evaluated.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">name2</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The name of the line being evaluated as specified by the genotyping lab. For lines that were not genotyped, this value is NA.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">studyName</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The experiment name.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">plotNumber</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The plot number within the field experiment</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">observationUnitName</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The name of the plot. This is the combination of experiment name and plot number.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">studyDesign</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The statistical design of the field experiment.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">replicate</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The replicate number within the experimental design.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">blockNumber</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The block number within the experimental design. This is eqivalent to the replicate number in this dataset.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">rowNumber</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The row position within the experiment where the experimental unit is located.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">colNumber</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;The column position within the experiment where the experimental unit is located.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">DON</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;Observation data for Deoxynivalenol content</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">FDK_V</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;Observation data for Fusarium Damaged Kernels (FDK) determined using conventional visual assessments</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">FDK_L</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;Observation data for FDK determined using manually labeled images.</span></span></span></p>
			</td>
		</tr>
		<tr>
			<td style="border-bottom:1px solid black; border-left:1px solid black; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:124px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">FDK_Lhat</span></span></span></p>
			</td>
			<td style="border-bottom:1px solid black; border-left:none; border-right:1px solid black; border-top:none; height:21px; vertical-align:bottom; width:612px">
			<p><span style="color:#000000"><span style="font-size:12px"><span style="font-family:Verdana,Geneva,sans-serif">&nbsp;Predicted FDK determined using predicted labels based on a neural network.</span></span></span></p>
			</td>
		</tr>
	</tbody>
</table>
