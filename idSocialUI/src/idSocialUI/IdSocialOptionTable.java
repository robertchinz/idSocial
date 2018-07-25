package idSocialUI;

import java.awt.EventQueue;

import net.miginfocom.swing.MigLayout;

import javax.swing.table.DefaultTableModel;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.EventObject;
import java.util.Hashtable;

import javax.swing.event.CellEditorListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;

import javax.swing.table.*;
import java.util.*;

import javax.swing.*;
import javax.swing.table.*;
import java.util.Vector;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;

public class IdSocialOptionTable extends JTableX {

	private JTableX table = new JTableX();
	private DefaultTableModel model;
	private RowEditorModel rm;

	public JTable getTrajectoryTable() {
		return table;
	}

	public IdSocialOptionTable() {
		super();
		rm = new RowEditorModel();
		this.setRowEditorModel(rm);
		model = new DefaultTableModel() {
			public boolean isCellEditable(int row, int column) {
				if (column == 1) {
					return true;
				}
				return false;
			}
		};
		this.setModel(model);

		// this.setDefaultRenderer(Double.class, new ColorRenderer(true));
		// table.setDefaultEditor(Value.class, new ValueEditor());
	}

	public IdSocialOptionTable(final Object[][] rowData, final Object[] colNames) {
		super(rowData, colNames);
		rm = new RowEditorModel();
		this.setRowEditorModel(rm);
		// rm = null;
	};

	public void setTableData(String[] parNames, Object[] parValues) {

		DefaultTableModel model = (DefaultTableModel) this.getModel();
		RowEditorModel rm = this.getRowEditorModel();

		model.setColumnIdentifiers(new Object[] { "Parameter", "Value" });
		// TableModel model = table.gemodelodel();
		// RowEditorModel rm = table.getRowEditorModel();

		for (int actRow = 0; actRow <= parNames.length - 1; actRow++) {

			if (model.getRowCount() <= actRow) {
				model.addRow(new Object[] { "", "" });
			}
			model.setValueAt(parNames[actRow], actRow, 0);

			Object elem = parValues[actRow];
			if (elem != null) {
				Class actClass = elem.getClass();
				if (elem instanceof String[]) {
					String[] combolist = (String[]) elem;
					model.setValueAt(combolist[0], actRow, 1);

					JComboBox cb = new JComboBox(combolist);
					// cb.setSelectedItem(1);
					((JLabel) cb.getRenderer()).setHorizontalAlignment(SwingConstants.CENTER);
					DefaultCellEditor ed = new DefaultCellEditor(cb);
					// inform the RowEditorMode of the situation
					rm.addEditorForRow(actRow, ed);
					// System.out.println(actRow);
					// System.out.println(elem.getClass());
					// System.out.println(ed);

				} else if (elem instanceof Boolean) {
					model.setValueAt(((Boolean) elem).booleanValue(), actRow, 1);
					JCheckBox chkBox = new JCheckBox("", ((Boolean) elem).booleanValue());
					chkBox.setHorizontalAlignment(JCheckBox.CENTER);
					DefaultCellEditor ed = new DefaultCellEditor(chkBox);
					// inform the RowEditorMode of the situation
					// System.out.println(actRow);
					// System.out.println(elem.getClass());
					// System.out.println(ed);

					rm.addEditorForRow(actRow, ed);

				} else if (elem instanceof double[]) {
					model.setValueAt(((double[]) elem), actRow, 1);

					// String doubleTxt = Arrays.toString((double[]) elem);
					// JTextField txtFld = new JTextField(doubleTxt);
					// DefaultCellEditor ed = new DefaultCellEditor(txtFld);
					
					DoubleArrayEditor ed = new DoubleArrayEditor();

					// inform the RowEditorMode of the situation
					System.out.println(actRow);
					System.out.println(elem.getClass());
					System.out.println(model.getValueAt(actRow, 1));
					// System.out.println(txtFld);
					// System.out.println(doubleTxt);

					 rm.addEditorForRow(actRow, ed);

				} else if (elem instanceof double[][]) {
					model.setValueAt((double[][]) elem, actRow, 1);
					this.setValueAt(Arrays.toString((double[][]) elem), actRow, 1);

					JTextField txtFld = new JTextField(Arrays.toString((double[][]) elem));
					// txtFld.setHorizontalAlignment(JCheckBox.CENTER);
					DefaultCellEditor ed = new DefaultCellEditor(txtFld);

					// inform the RowEditorMode of the situation
					// System.out.println(actRow);
					// System.out.println(elem.getClass());
					// System.out.println(ed);
					// System.out.println(Arrays.toString((double[][]) elem));

				} else {
					model.setValueAt(parValues[actRow], actRow, 1);
					System.out.println(elem.getClass());

				}

			}

		}
		this.setModel(model); // table.setModel(model);

	}

	public Object[] getTableData() {

		int rows = this.getRowCount();
		DefaultTableModel model = (DefaultTableModel) this.getModel();

		Object[] colData = new Object[rows];

		if (rows > 0) {
			for (int actRow = 0; actRow <= rows - 1; actRow++) {
				colData[actRow] = model.getValueAt(actRow, 1);
				System.out.println(colData[actRow]);
			}
		}
		return colData;
	}

}

class DoubleArrayRenderer extends DefaultTableCellRenderer {

	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,
			int row, int column) {

		value = Arrays.toString((double[]) value);

		// And pass it on to parent class

		return super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
	}
}

class DoubleArrayEditor extends DefaultCellEditor {

public DoubleArrayEditor() {
    super(new JTextField());
  }

  public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected,
      int row, int column) {
    JTextField editor = (JTextField) super.getTableCellEditorComponent(table, Arrays.toString((double[]) value), isSelected,
        row, column);

    
    return editor;
  }
}
  
// http://www.javaworld.com/article/2077465/learn-java/java-tip-102--add-multiple-jtable-cell-editors-per-column.html
class RowEditorModel {
	private Hashtable data;

	public RowEditorModel() {
		data = new Hashtable();
	}

	public void addEditorForRow(int row, TableCellEditor e) {
		data.put(new Integer(row), e);
	}

	public void removeEditorForRow(int row) {
		data.remove(new Integer(row));
	}

	public TableCellEditor getEditor(int row) {
		return (TableCellEditor) data.get(new Integer(row));
	}
}

class JTableX extends JTable {
	protected RowEditorModel rm;

	public JTableX() {
		super();
		rm = null;
	}

	public JTableX(TableModel tm) {
		super(tm);
		rm = null;
	}

	public JTableX(TableModel tm, TableColumnModel cm) {
		super(tm, cm);
		rm = null;
	}

	public JTableX(TableModel tm, TableColumnModel cm, ListSelectionModel sm) {
		super(tm, cm, sm);
		rm = null;
	}

	public JTableX(int rows, int cols) {
		super(rows, cols);
		rm = null;
	}

	public JTableX(final Vector rowData, final Vector columnNames) {
		super(rowData, columnNames);
		rm = null;
	}

	public JTableX(final Object[][] rowData, final Object[] colNames) {
		super(rowData, colNames);
		rm = null;
	}

	// new constructor
	public JTableX(TableModel tm, RowEditorModel rm) {
		super(tm, null, null);
		this.rm = rm;
	}

	public void setRowEditorModel(RowEditorModel rm) {
		this.rm = rm;
	}

	public RowEditorModel getRowEditorModel() {
		return rm;
	}

	public TableCellEditor getCellEditor(int row, int col) {
		TableCellEditor tmpEditor = null;
		if (rm != null)
			tmpEditor = rm.getEditor(row);
		if (tmpEditor != null)
			return tmpEditor;
		return super.getCellEditor(row, col);

	}

	public TableCellRenderer getCellRenderer(int row, int column) {
		if (column == 1) {
			Object value = getValueAt(row, column);
			if (value instanceof Boolean) {
				return getDefaultRenderer(Boolean.class);
			}
			// The following was done by me
			// if (value instanceof Double) {
			// return getDefaultRenderer(Double.class);
			// }
			if (value instanceof double[]) {
				// return getDefaultRenderer(Double[].class);
				TableCellRenderer DArrayRenderer = new DoubleArrayRenderer();
				return DArrayRenderer;
				// ColorRenderer cr = new ColorRenderer(true);
				// return cr;
			}
			// end
		}
		return super.getCellRenderer(row, column);
	}
}
