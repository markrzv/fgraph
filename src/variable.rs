use std::cell::RefCell;
use std::rc::Rc;

use crate::gaussian::Gaussian;

#[derive(Clone)]
pub struct Variable {
    value: Rc<RefCell<Gaussian>>,
}

impl Variable {
    pub fn new(value: Gaussian) -> Self {
        Self { value: Rc::new(RefCell::new(value)) }
    }

    pub fn default() -> Self {
        Self::new(Gaussian::default())
    }

    pub fn get(&self) -> Gaussian {
        *self.value.borrow()
    }

    pub fn set(&self, new_value: Gaussian) {
        *self.value.borrow_mut() = new_value;
    }
}