use noodles::bam::Record;

pub trait AddRecord {
    fn add_record(&mut self, record: &Record);
}