use std::vec::Vec;

#[derive(Copy, Clone)]
pub struct Interval<T: Sized+Copy> {
    pub start: i32,
    pub stop: i32,
    pub value: T
}

impl<T: Sized+Copy> Interval<T> {
    pub fn new(s: i32, e: i32, v: T) -> Interval<T> {
        Interval::<T>{
            start: s,
            stop: e,
            value: v
        }
    }
}

#[derive(Copy, Clone)]
pub struct IntervalTreeOpts {
    pub depth: u32,
    pub minbucket: usize,
    pub leftextent: i32,
    pub rightextent: i32,
    pub maxbucket: usize
}

pub struct IntervalTree<T: Sized+Copy> {
    intervals: Vec<Interval<T>>,
    left: Option<Box<IntervalTree<T>>>,
    right: Option<Box<IntervalTree<T>>>,
    center: i32
}

impl<T: Sized+Copy> IntervalTree<T> {
    pub fn new() -> IntervalTree<T> {
        IntervalTree::<T>{
            intervals: Vec::new(),
            left: None,
            right: None,
            center: 0
        }
    }

    pub fn new_from(ivals: &Vec<Interval<T>>) -> IntervalTree<T> {
        Self::new_from_opts(ivals, &IntervalTreeOpts{
            depth: 16,
            minbucket: 64,
            leftextent: 0,
            rightextent: 0,
            maxbucket: 512
        })
    }
    pub fn new_from_opts(ivals: &Vec<Interval<T>>, opts: &IntervalTreeOpts) -> IntervalTree<T> {
        let mut this = Self::new();
        let ref mut ivals = ivals.clone();
        let ref mut opts = opts.clone();
        opts.depth -= 1;
        if opts.depth == 0 || (ivals.len() < opts.minbucket && ivals.len() < opts.maxbucket) {
            this.intervals = ivals.clone();
        }
        else {
            if opts.leftextent == 0 && opts.rightextent == 0 {
                ivals.sort_by(|a, b| a.start.cmp(&b.start));
            }

            let leftp: i32;
            let mut rightp: i32 = 0;
            let centerp: i32;

            if opts.leftextent != 0 || opts.rightextent != 0 {
                leftp = opts.leftextent;
                rightp = opts.rightextent;
            }
            else {
                leftp = ivals[0].start;
                let mut stops: Vec<i32> = Vec::new();
                stops.resize(ivals.len(), 0);
                for i in 0..ivals.len() {
                    stops[i] = ivals[i].stop;
                    if rightp < stops[i] {rightp = stops[i]}
                }
            }

            // centerp = (leftp + rightp) / 2;
            centerp = ivals[ivals.len() / 2].start;
            this.center = centerp;

            let mut lefts: Vec<Interval<T>> = Vec::new();
            let mut rights: Vec<Interval<T>> = Vec::new();

            for i in 0..ivals.len() {
                let interval: Interval<T> = ivals[i];
                if interval.stop < this.center {
                    lefts.push(interval);
                }
                else if interval.start > this.center {
                    rights.push(interval);
                }
                else {
                    this.intervals.push(interval);
                }
            }

            if !lefts.is_empty() {
                this.left = Some(Box::new(IntervalTree::new_from_opts(&mut lefts, &IntervalTreeOpts{
                    depth: opts.depth,
                    minbucket: opts.minbucket,
                    leftextent: leftp,
                    rightextent: centerp,
                    maxbucket: opts.maxbucket
                })));
            }
            if !rights.is_empty() {
                this.right = Some(Box::new(IntervalTree::new_from_opts(&mut rights, &IntervalTreeOpts{
                    depth: opts.depth,
                    minbucket: opts.minbucket,
                    leftextent: centerp,
                    rightextent: rightp,
                    maxbucket: opts.maxbucket
                })));
            }
        }
        this
    }

    pub fn find_overlapping(&self, start: i32, stop: i32, overlapping: &mut Vec<Interval<T>>) {
        if !self.intervals.is_empty() && !(stop < self.intervals[0].start) {
            for i in 0..self.intervals.len() {
                let interval = self.intervals[i];
                if interval.stop >= start && interval.start <= stop {
                    overlapping.push(interval);
                }
            }
        }

        if self.left.is_some() && start <= self.center  {
            self.left.as_ref().unwrap().find_overlapping(start, stop, overlapping);
        }
        if self.right.is_some() && stop >= self.center {
            self.right.as_ref().unwrap().find_overlapping(start, stop, overlapping);
        }
    }

    #[allow(dead_code)]
    pub fn find_contained(&self, start: i32, stop: i32, contained: &mut Vec<Interval<T>>) {
        if !self.intervals.is_empty() && !(stop < self.intervals[0].start) {
            for i in 0..self.intervals.len() {
                let interval = self.intervals[i];
                if interval.start >= start && interval.stop <= stop {
                    contained.push(interval);
                }
            }
        }

        if self.left.is_some() && start <= self.center  {
            self.left.as_ref().unwrap().find_contained(start, stop, contained);
        }
        if self.right.is_some() && stop >= self.center {
            self.right.as_ref().unwrap().find_contained(start, stop, contained);
        }
    }
}
impl<T: Sized+Copy> Clone for IntervalTree<T> {
    fn clone(&self) -> Self {
        IntervalTree::<T> {
            intervals: self.intervals.clone(),
            left: self.left.clone(),
            right: self.right.clone(),
            center: self.center
        }
    }

    fn clone_from(&mut self, source: &Self) {
        self.center = source.center;
        self.intervals = source.intervals.clone();
        if source.left.is_some() {
            self.left = source.left.clone();
        }
        else {
            self.left = None;
        }
        if source.right.is_some() {
            self.right = source.right.clone();
        }
        else {
            self.right = None;
        }
    }
}
