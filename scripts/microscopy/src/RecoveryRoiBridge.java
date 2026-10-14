import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.io.*;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.process.*;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.file.*;
import java.util.*;
import java.util.zip.*;
import javax.imageio.ImageIO;

/** Decode and rasterize bounded ROI exchanges using pinned ImageJ semantics. */
public final class RecoveryRoiBridge {
    static byte[] raster(Roi roi, int w, int h) {
        int type = roi.getType();
        if (type != Roi.RECTANGLE && type != Roi.POLYGON && type != Roi.FREEROI && type != Roi.COMPOSITE)
            throw new IllegalArgumentException("Untested ROI type: " + roi.getTypeAsString());
        if ((type == Roi.RECTANGLE && (roi.getClass() != Roi.class || roi.subPixelResolution() || roi.getCornerDiameter() != 0))
                || ((type == Roi.POLYGON || type == Roi.FREEROI) && (roi.getClass() != PolygonRoi.class || ((PolygonRoi) roi).isSplineFit()))
                || (type == Roi.COMPOSITE && roi.getClass() != ShapeRoi.class))
            throw new IllegalArgumentException("Untested ROI subtype, rounded/subpixel rectangle, or spline");
        if (!roi.isArea() || roi.getPosition() != 0 || roi.getCPosition() != 0 || roi.getZPosition() != 0 || roi.getTPosition() != 0)
            throw new IllegalArgumentException("Unpositioned area ROI required");
        Rectangle b = roi.getBounds();
        if (b.x < 0 || b.y < 0 || b.width <= 0 || b.height <= 0 || (long)b.x + b.width > w || (long)b.y + b.height > h)
            throw new IllegalArgumentException("ROI bounds outside crop");
        ImageProcessor local = roi.getMask();
        if (local == null && type != Roi.RECTANGLE)
            throw new IllegalArgumentException("Non-rectangle lacks ImageJ mask");
        if (local != null && (local.getWidth() != b.width || local.getHeight() != b.height))
            throw new IllegalArgumentException("Local mask dimensions differ from bounds");
        byte[] mask = new byte[Math.multiplyExact(w, h)];
        int count = 0;
        for (int y = b.y; y < b.y + b.height; y++)
            for (int x = b.x; x < b.x + b.width; x++)
                if (local == null || local.get(x-b.x, y-b.y) != 0) { mask[y*w+x] = (byte)255; count++; }
        if (count == 0) throw new IllegalArgumentException("Empty ROI mask");
        return mask;
    }

    static void png(Path path, byte[] mask, int w, int h) throws Exception {
        if (Files.exists(path, LinkOption.NOFOLLOW_LINKS)) throw new IOException("Existing output: " + path);
        BufferedImage image = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY);
        image.getRaster().setDataElements(0, 0, w, h, mask);
        if (!ImageIO.write(image, "PNG", path.toFile())) throw new IOException("PNG writer missing");
    }

    static Roi decodeBytes(byte[] encoded, String name) throws Exception {
        if (encoded.length < 64 || encoded[0] != 'I' || encoded[1] != 'o' || encoded[2] != 'u' || encoded[3] != 't'
                || ByteBuffer.wrap(encoded).getShort(RoiDecoder.SUBTYPE) != 0)
            throw new IllegalArgumentException("Malformed ROI header or untested encoded subtype");
        Roi roi = new RoiDecoder(encoded, name).getRoi();
        if (roi == null) throw new IllegalArgumentException("Decoder returned null");
        // ZIP entry identity owns the role; encoded display names never redirect it.
        return roi;
    }

    static void decode(Path archive, Path out, int w, int h) throws Exception {
        if (w <= 0 || h <= 0 || (long)w*h > 100000000) throw new IllegalArgumentException("Bounded dimensions required");
        if (!Files.isDirectory(out) || Files.isSymbolicLink(out)) throw new IOException("Fresh caller staging required");
        List<String> records = new ArrayList<String>();
        Set<String> names = new HashSet<String>();
        int total = 0;
        try (ZipFile zip = new ZipFile(archive.toFile())) {
            Enumeration<? extends ZipEntry> entries = zip.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement(); String name = entry.getName();
                if (entry.isDirectory() || (!name.equals("nucleus.roi") && !name.matches("background_[0-9]{3}\\.roi"))
                        || name.equals("background_000.roi") || !names.add(name) || names.size() > 64)
                    throw new IllegalArgumentException("Invalid or duplicate ROI role: " + name);
                ByteArrayOutputStream bytes = new ByteArrayOutputStream();
                try (InputStream in = zip.getInputStream(entry)) {
                    byte[] block = new byte[8192]; int n;
                    while ((n = in.read(block)) != -1) {
                        bytes.write(block, 0, n); total += n;
                        if (bytes.size() > 1048576 || total > 16777216) throw new IllegalArgumentException("Oversized ROI archive");
                    }
                }
                Roi roi = decodeBytes(bytes.toByteArray(), name);
                byte[] mask = raster(roi, w, h);
                String stem = name.substring(0, name.length()-4);
                png(out.resolve(stem + ".png"), mask, w, h);
                long count = 0; for (byte b : mask) if (b != 0) count++;
                Rectangle bounds = roi.getBounds();
                records.add(stem + "\t" + roi.getTypeAsString() + "\t" + count + "\t" + bounds.x + "\t" + bounds.y + "\t" + bounds.width + "\t" + bounds.height);
            }
        }
        if (!names.contains("nucleus.roi")) throw new IllegalArgumentException("Exactly one nucleus.roi required");
        Files.write(out.resolve("masks.tsv"), records, java.nio.charset.StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
    }

    static void zip(Path path, Roi roi) throws Exception {
        try (ZipOutputStream zip = new ZipOutputStream(Files.newOutputStream(path, StandardOpenOption.CREATE_NEW))) {
            zip.putNextEntry(new ZipEntry("nucleus.roi")); zip.write(RoiEncoder.saveAsByteArray(roi)); zip.closeEntry();
        }
    }

    static void fixture(Path out) throws Exception {
        Files.createDirectory(out);
        int w = 16, h = 12;
        short[] pixels = new short[w*h];
        for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) pixels[y*w+x] = (short)(13+x*237+y*1051);
        short[] before = pixels.clone();
        ImagePlus original = new ImagePlus("synthetic", new ShortProcessor(w,h,pixels,null));
        if (!new FileSaver(original).saveAsTiff(out.resolve("original_before.tif").toString())) throw new IOException("TIFF write failed");
        LinkedHashMap<String,Roi> cases = new LinkedHashMap<String,Roi>();
        cases.put("asymmetric_integer",new PolygonRoi(new int[]{1,10,12,6,2},new int[]{1,2,7,10,7},5,Roi.POLYGON));
        cases.put("concave_freehand",new PolygonRoi(new int[]{2,13,13,8,8,5,5,2},new int[]{2,2,10,10,5,5,10,10},8,Roi.FREEROI));
        cases.put("subpixel_polygon",new PolygonRoi(new float[]{1.2f,10.7f,12.1f,6.25f,2.8f},new float[]{1.8f,2.2f,7.6f,10.3f,7.15f},5,Roi.POLYGON));
        cases.put("composite_hole",new ShapeRoi(new Roi(1,1,13,10)).not(new ShapeRoi(new Roi(5,4,4,3))));
        cases.put("boundary_rectangle",new Roi(0,0,4,h));
        List<String> records = new ArrayList<String>();
        for (Map.Entry<String,Roi> entry : cases.entrySet()) {
            String label = entry.getKey(); Roi roi = entry.getValue(); roi.setName("encoded_display_name");
            zip(out.resolve(label+".zip"),roi);
            Roi decoded = decodeBytes(RoiEncoder.saveAsByteArray(roi), "nucleus.roi");
            byte[] mask = raster(decoded,w,h); long sum = 0, count = 0;
            for (int i = 0; i < mask.length; i++) if (mask[i] != 0) { count++; sum += original.getProcessor().get(i); }
            original.setRoi(decoded); ResultsTable results = new ResultsTable();
            new Analyzer(original,Measurements.INTEGRATED_DENSITY,results).measure();
            double raw = results.getValue("RawIntDen",0);
            if (raw != sum) throw new AssertionError("ImageJ RawIntDen differs from selected pixel sum");
            png(out.resolve(label+".png"),mask,w,h);
            records.add(label+"\t"+count+"\t"+sum+"\t"+raw);
        }
        original.deleteRoi();
        if (!Arrays.equals(before,(short[])original.getProcessor().getPixels())) throw new AssertionError("Source pixels changed");
        if (!new FileSaver(original).saveAsTiff(out.resolve("original_after.tif").toString())) throw new IOException("TIFF write failed");
        cases.clear(); cases.put("line",new Line(1,1,8,7)); cases.put("point",new PointRoi(3,4));
        Roi positioned = new Roi(2,2,4,4); positioned.setPosition(0,2,0); cases.put("z",positioned);
        positioned = new Roi(2,2,4,4); positioned.setPosition(2); cases.put("stack",positioned);
        cases.put("rotated",new RotatedRectRoi(2,2,8,7,3));
        Roi rounded = new Roi(2,2,8,7); rounded.setCornerDiameter(3); cases.put("rounded",rounded);
        PolygonRoi spline = new PolygonRoi(new int[]{1,10,12,6,2},new int[]{1,2,7,10,7},5,Roi.POLYGON); spline.fitSpline(); cases.put("spline",spline);
        cases.put("out_of_bounds",new Roi(-1,1,4,4)); cases.put("subpixel_rectangle",new Roi(1.2,1.1,3.3,4.4));
        cases.put("ellipse",new OvalRoi(1,1,4,4));
        for (Map.Entry<String,Roi> entry : cases.entrySet()) zip(out.resolve("reject_"+entry.getKey()+".zip"),entry.getValue());
        Files.write(out.resolve("fixture.tsv"),records,java.nio.charset.StandardCharsets.UTF_8,StandardOpenOption.CREATE_NEW);
        System.out.println("Five genuine ImageJ geometries: RawIntDen equality and unchanged source pixels passed");
    }

    public static void main(String[] args) throws Exception {
        System.out.println("ImageJ="+IJ.getVersion()+" Java="+System.getProperty("java.version")+" headless="+System.getProperty("java.awt.headless"));
        if (!IJ.getVersion().equals("1.54p")) throw new IllegalStateException("ImageJ version mismatch");
        if (args.length == 2 && args[0].equals("fixture")) fixture(Paths.get(args[1]));
        else if (args.length == 5 && args[0].equals("raster")) decode(Paths.get(args[1]),Paths.get(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]));
        else throw new IllegalArgumentException("fixture OUT or raster ZIP OUT WIDTH HEIGHT");
    }
}
