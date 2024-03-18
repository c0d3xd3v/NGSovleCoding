import QtQuick 6.6
import QtQuick.Controls 6.6
import QtQuick.Layouts 6.6

ColumnLayout {
    id: columnLayout
    anchors.fill: parent
    spacing: 10

    Label {
        id: label
        text: qsTr("Meshing")
        horizontalAlignment: Text.AlignHCenter
        verticalAlignment: Text.AlignVCenter
        wrapMode: Text.WordWrap
        renderTypeQuality: Text.HighRenderTypeQuality
        font.underline: false
        font.italic: false
        font.bold: false
        Layout.fillWidth: true
    }

    RowLayout {
        id: rowLayout
        width: 100
        height: 100

        HidableLable {}


        /*
        Label {
            id: label2
            text: qsTr("Label")
        }
        */
        TextField {
            id: textField
            Layout.fillWidth: true
            placeholderText: qsTr("Text Field")
        }
    }

    RowLayout {
        id: rowLayout1
        width: 100
        height: 100

        HidableLable {}

        TextField {
            id: textField1
            Layout.fillWidth: true
            placeholderText: qsTr("Text Field")
        }
    }

    RowLayout {
        id: rowLayout2
        width: 100
        height: 100

        HidableLable {}

        TextField {
            id: textField2
            Layout.fillWidth: true
            placeholderText: qsTr("Text Field")
        }
    }

    Button {
        id: button1
        text: qsTr("Button")
        Layout.fillWidth: true
    }

    ToolSeparator {
        id: toolSeparator
        Layout.fillWidth: true
        orientation: Qt.Horizontal
    }

    CheckBox {
        id: checkBox
        text: qsTr("Check Box")
        Layout.fillWidth: true
    }

    Slider {
        id: slider
        value: 0.5
        Layout.fillWidth: true
    }

    Switch {
        id: switch1
        text: qsTr("Switch")
        Layout.fillWidth: true
    }

    Item {
        id: item1
        width: 200
        height: 200
        Layout.fillHeight: true
        Layout.fillWidth: true
    }
}
