import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

// One builder; previous exports and builder are preserved in this run.
const run = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const out = path.join(run, 'tables/integrated');
const qa = path.join(run, 'qa/integrated_tables');
await fs.mkdir(qa, {recursive:true});
const payload = JSON.parse(await fs.readFile(path.join(out, 'workbook_payload.json'), 'utf8'));
const workbook = Workbook.create();
const font = {name:'Arial', size:10, color:'#222222'};
function col(index) {
  let result = '';
  for (let n=index+1; n>0; n=Math.floor((n-1)/26)) result=String.fromCharCode(65+(n-1)%26)+result;
  return result;
}
for (const spec of payload.sheets) {
  const sheet = workbook.worksheets.add(spec.name);
  const lastCol = col(spec.headers.length-1);
  const lastRow = 5 + spec.rows.length;
  sheet.showGridLines = false;
  sheet.getRange(`A1:${lastCol}${lastRow}`).format = {font,verticalAlignment:'center',rowHeight:25};
  sheet.getRange('A2').values = [[spec.title]];
  sheet.getRange('A2').format.font = {...font,size:14,bold:true};
  sheet.getRange(`A3:${lastCol}3`).format.borders = {bottom:{style:'thin',color:'#7E8E99'}};
  sheet.getRange(`A5:${lastCol}5`).values = [spec.headers];
  sheet.getRange(`A5:${lastCol}5`).format = {font:{...font,bold:true},fill:'#E8EDF0',wrapText:true,
    horizontalAlignment:'center',verticalAlignment:'center',rowHeight:42,
    borders:{bottom:{style:'thin',color:'#7E8E99'}}};
  sheet.getRange(`A6:${lastCol}${lastRow}`).values = spec.rows.map(row=>row.map(v=>v===null?'n.a.  ':v));
  sheet.getRange(`A6:${lastCol}${lastRow}`).format.wrapText = true;
  sheet.getRange(`A6:${lastCol}${lastRow}`).format.rowHeight = spec.name==='Species measurements'?32:45;
  for (let c=0;c<spec.headers.length;c++) {
    sheet.getRange(`${col(c)}1:${col(c)}${lastRow}`).format.columnWidth = spec.widths[c];
    if (String(c) in spec.formats) {
      const range = sheet.getRange(`${col(c)}6:${col(c)}${lastRow}`);
      range.setNumberFormat(spec.formats[String(c)] + '"  "');
      range.format.horizontalAlignment = 'right';
    }
  }
  sheet.getRange(`A${lastRow}:${lastCol}${lastRow}`).format.borders = {bottom:{style:'thin',color:'#7E8E99'}};
  // Definitions sit beside the table without altering its data or row heights.
  const noteCol=col(spec.headers.length+1);
  sheet.getRange(`${noteCol}5`).values=[['Definitions and interpretation']];
  sheet.getRange(`${noteCol}5`).format.font={...font,bold:true};
  const sentences=spec.caption.replaceAll('D. fuscus','D~fuscus').split(/(?<=[.!?])\s+/).map(s=>s.replaceAll('D~fuscus','D. fuscus'));
  const noteRows=sentences.map(s=>[s.trim()]);
  sheet.getRange(`${noteCol}5:${noteCol}${5+noteRows.length}`).format.columnWidth=100;
  sheet.getRange(`${noteCol}6:${noteCol}${5+noteRows.length}`).values=noteRows;
  sheet.getRange(`${noteCol}6:${noteCol}${5+noteRows.length}`).format={font:{...font,color:'#555555'},wrapText:true,verticalAlignment:'center'};
  if (spec.rows.length>15) {
    sheet.freezePanes.freezeRows(5);
    if (spec.headers.length>9) sheet.freezePanes.freezeColumns(1);
  }
  sheet.tabColor='#0072B2';
}
workbook.recalculate();
const inspections=[];
for (const spec of payload.sheets) {
  const lastCol=col(spec.headers.length-1);
  const check=await workbook.inspect({kind:'table',range:`'${spec.name}'!A5:${lastCol}${Math.min(10,5+spec.rows.length)}`,
    include:'values,formulas',tableMaxRows:6,tableMaxCols:spec.headers.length,maxChars:7000});
  inspections.push(check.ndjson);
  const preview=await workbook.render({sheetName:spec.name,range:`A1:${lastCol}${Math.min(14,5+spec.rows.length)}`,scale:1.2,format:'png'});
  await fs.writeFile(path.join(qa,`${spec.name.replaceAll(' ','_')}.png`),new Uint8Array(await preview.arrayBuffer()));
  const noteCol=col(spec.headers.length+1);
  const notesPreview=await workbook.render({sheetName:spec.name,range:`${noteCol}5:${noteCol}12`,scale:1.2,format:'png'});
  await fs.writeFile(path.join(qa,`${spec.name.replaceAll(' ','_')}_definitions.png`),new Uint8Array(await notesPreview.arrayBuffer()));
  if (spec.rows.length>9) {
    const end=5+spec.rows.length;
    const endPreview=await workbook.render({sheetName:spec.name,range:`A${Math.max(15,end-8)}:${lastCol}${end}`,scale:1.2,format:'png'});
    await fs.writeFile(path.join(qa,`${spec.name.replaceAll(' ','_')}_end.png`),new Uint8Array(await endPreview.arrayBuffer()));
  }
}
await fs.writeFile(path.join(qa,'inspection.ndjson'),inspections.join('\n'));
const errors=await workbook.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',
  options:{useRegex:true,maxResults:100},summary:'Final formula error scan'});
await fs.writeFile(path.join(qa,'error_scan.ndjson'),errors.ndjson);
const xlsx=await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(path.join(out,'Manuscript_tables_integrated.xlsx'));
console.log(`Created ${payload.sheets.length} consolidated tables; original workbook retained.`);
