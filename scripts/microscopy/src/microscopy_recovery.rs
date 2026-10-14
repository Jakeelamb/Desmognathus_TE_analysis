use std::io::{self, Read};

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 { return Err("usage: microscopy_recovery BACKGROUND < uint8_pixels".into()); }
    let background: f64 = args[1].parse().map_err(|_| "invalid background")?;
    if !background.is_finite() || !(1.0..=255.0).contains(&background) {
        return Err("background must be finite in [1,255]".into());
    }
    let mut pixels = Vec::new();
    io::stdin().take(100_000_001).read_to_end(&mut pixels).map_err(|e| e.to_string())?;
    if pixels.is_empty() || pixels.len() > 100_000_000 {
        return Err("nonempty bounded uint8 vector required".into());
    }
    let mut raw_sum: u64 = 0;
    let (mut min, mut max) = (255u8, 0u8);
    let (mut zeros, mut ones, mut saturated) = (0u64, 0u64, 0u64);
    let mut histogram = [0u64; 256];
    let (mut sum, mut correction) = (0.0f64, 0.0f64);
    for &pixel in &pixels {
        raw_sum += pixel as u64;
        min = min.min(pixel); max = max.max(pixel);
        zeros += u64::from(pixel == 0); ones += u64::from(pixel == 1);
        saturated += u64::from(pixel == 255);
        histogram[usize::from(pixel)] += 1;
        let term = (background / f64::from(pixel.max(1))).log10();
        let next = sum + term;
        correction += if sum.abs() >= term.abs() { (sum - next) + term } else { (term - next) + sum };
        sum = next;
    }
    let iod = sum + correction;
    let lower_rank = (pixels.len() as u64 - 1) / 2;
    let upper_rank = pixels.len() as u64 / 2;
    let (mut cumulative, mut lower, mut upper) = (0u64, 0u8, 0u8);
    for (pixel, count) in histogram.iter().enumerate() {
        let next = cumulative + count;
        if cumulative <= lower_rank && next > lower_rank { lower = pixel as u8; }
        if cumulative <= upper_rank && next > upper_rank { upper = pixel as u8; break; }
        cumulative = next;
    }
    let median = (f64::from(lower) + f64::from(upper)) / 2.0;
    println!("{{\"count\":{},\"raw_sum\":{},\"min_pixel\":{},\"max_pixel\":{},\"zero_pixels\":{},\"one_pixels\":{},\"saturated_pixels\":{},\"background\":{},\"iod\":{},\"mean_od\":{},\"median_pixel\":{}}}",
        pixels.len(), raw_sum, min, max, zeros, ones, saturated, background, iod, iod / pixels.len() as f64, median);
    Ok(())
}

fn main() {
    if let Err(message) = run() { eprintln!("{}", message); std::process::exit(2); }
}
