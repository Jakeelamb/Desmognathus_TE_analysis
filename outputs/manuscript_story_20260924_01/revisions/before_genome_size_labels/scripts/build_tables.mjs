import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const run = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const payload = JSON.parse(await fs.readFile(path.join(run, 'tables/tables.json'), 'utf8'));
const workbook = Workbook.create();
const coverage = workbook.worksheets.add('Data coverage');
const associations = workbook.worksheets.add('Associations');

function base(sheet, lastCol, endRow) {
  sheet.showGridLines = false;
  sheet.getRange(`A1:${lastCol}${endRow}`).format.font = {name: 'Arial', size: 10, color: '#222222'};
  sheet.getRange(`A1:${lastCol}${endRow}`).format.verticalAlignment = 'center';
  sheet.getRange(`A5:${lastCol}5`).format = {fill: '#E8EDF0', font: {name:'Arial',size:10,bold:true,color:'#222222'}, wrapText:true, rowHeight:30};
  sheet.getRange(`A5:${lastCol}5`).format.borders = {bottom:{style:'thin',color:'#6F6F6F'}};
  sheet.getRange('A2').format.font = {name:'Arial',size:14,bold:true,color:'#222222'};
  sheet.getRange(`A3:${lastCol}3`).format.borders = {bottom:{style:'thin',color:'#6F6F6F'}};
  sheet.tabColor = '#0072B2';
}

base(coverage, 'F', 17);
coverage.getRange('A2').values = [['Table 1. Data coverage']];
coverage.getRange('A5:F5').values = [['Component','Species','Objects','Animals','Measurement','Source']];
coverage.getRange('A6:F11').values = payload.coverage.map(row => row.map(v => v === null ? 'n.a.' : v));
coverage.getRange('A6:F11').format.wrapText = true;
coverage.getRange('A6:F11').format.rowHeight = 44;
coverage.getRange('B6:D11').setNumberFormat('#,##0');
coverage.getRange('B6:D11').format.horizontalAlignment = 'center';
for (const [col,width] of Object.entries({A:38,B:10,C:12,D:11,E:43,F:25})) coverage.getRange(`${col}1:${col}17`).format.columnWidth = width;
coverage.mergeCells('A13:F16');
coverage.getRange('A13').values = [[payload.coverage_notes]];
coverage.getRange('A13:F16').format.wrapText = true;
coverage.getRange('A13:F16').format.verticalAlignment = 'top';
coverage.getRange('A13:F16').format.font = {name:'Arial',size:10,color:'#555555'};

base(associations, 'H', 22);
associations.getRange('A2').values = [['Table 2. Comparative associations']];
associations.getRange('A5:H5').values = [['Response or comparison','Species','Statistic','Estimate','95% CI low','95% CI high','P','Lambda']];
const data = payload.associations.map(r => [r.comparison,r.n_species,r.effect_type,r.estimate,r.ci_low ?? 'n.a.',r.ci_high ?? 'n.a.',r.p_value,r.lambda ?? 'n.a.']);
const lastRow = 5 + data.length;
associations.getRange(`A6:H${lastRow}`).values = data;
associations.getRange(`A6:H${lastRow}`).format.wrapText = true;
associations.getRange(`A6:H${lastRow}`).format.rowHeight = 46;
associations.getRange(`B6:B${lastRow}`).setNumberFormat('0');
associations.getRange(`D6:F${lastRow}`).setNumberFormat('0.0000');
associations.getRange(`G6:G${lastRow}`).setNumberFormat('0.00E+00');
associations.getRange(`H6:H${lastRow}`).setNumberFormat('0.00E+00');
associations.getRange(`B6:H${lastRow}`).format.horizontalAlignment = 'right';
associations.getRange(`B6:C${lastRow}`).format.horizontalAlignment = 'center';
for (const [col,width] of Object.entries({A:52,B:9,C:23,D:12,E:13,F:13,G:14,H:10})) associations.getRange(`${col}1:${col}22`).format.columnWidth = width;
for (let i=0; i<data.length; i++) {
  if (i === 0 || payload.associations[i].block !== payload.associations[i-1].block) associations.getRange(`A${i+6}:H${i+6}`).format.borders = {top:{style:'thin',color:'#9AA7AE'}};
}
associations.mergeCells('A16:H19');
associations.getRange('A16').values = [[payload.association_notes]];
associations.getRange('A16:H19').format.wrapText = true;
associations.getRange('A16:H19').format.verticalAlignment = 'top';
associations.getRange('A16:H19').format.font = {name:'Arial',size:10,color:'#555555'};
associations.mergeCells('A21:H22');
associations.getRange('A21').values = [['Sources: composition test export; Figure 3 full-sample model export; S26 cellular PGLS; retained LTR PGLS export. Full source identifiers and coefficient scales accompany the downloadable CSV.']];
associations.getRange('A21:H22').format.wrapText = true;
associations.getRange('A21:H22').format.font = {name:'Arial',size:10,color:'#555555'};

workbook.recalculate();
const inspect = await workbook.inspect({kind:'table',range:'Associations!A5:H13',include:'values,formulas',tableMaxRows:10,tableMaxCols:8,maxChars:6500});
await fs.writeFile(path.join(run,'qa/workbook_inspection.ndjson'),inspect.ndjson);
const errors = await workbook.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Final formula error scan'});
await fs.writeFile(path.join(run,'qa/workbook_error_scan.ndjson'),errors.ndjson);
for (const [sheetName,range,file] of [['Data coverage','A1:F17','Table_1_preview.png'],['Associations','A1:H22','Table_2_preview.png']]) {
  const preview = await workbook.render({sheetName,range,scale:1.5,format:'png'});
  await fs.writeFile(path.join(run,'qa',file),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(path.join(run,'tables/Manuscript_tables.xlsx'));
console.log('Created two editable manuscript tables; previews and inspection saved.');
